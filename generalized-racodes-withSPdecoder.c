#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define Q 3 //number of repetition (default 3)
#define A 1 //number of edges at combiner (this must be one for this program)
#define INFO_SIZE 5000 //number of information bits
#define CODE_SIZE (INFO_SIZE*Q/A) //must be integer
#define DETAIL    0 //show details of tanner graph
#define NONSYS    1 //nonsystematic mode (NOTE: A must be 1 for this option)

struct member{
  int    value;
  struct member *next;
};
void info_gen(int *data);
void random_interleaver(int *order);
void sys_ra_encoder(int *data, int *code, int *order);
void printintseq(int *a,int size);
double bpsk(int b);
double gauss(double av,double st);
void generate_tanner_graph(int *order,int vn[][Q],int cn[][A+2],int *vd,int *cd);
double llr_awgn(double y, double sig2);
void sum_product_decoder(double *llr,int *ans,int vn[][Q],int cn[][A+2],int *vd,int *cd, int itemax);
void printdoubleseq(double *a,int size);
void list_add(int tmp, struct member *list);
void list_del(int tmp, struct member *list);
void list_gen(int n, struct member *list);
int  list_display(int i, struct member *list);
void list_all_display(struct member *list);
void list_release(struct member *list);
void ra_random_interleaver(int *order);
void show_graph_details(int *data, int *code, int *order, int *vd, int *cd, int vn[][Q],int cn[][A+2]);

int main(int argc, char *argv[]){
  int data[INFO_SIZE],code[CODE_SIZE+INFO_SIZE],order[INFO_SIZE*Q],i,ans[INFO_SIZE],itemax;
  int vn[INFO_SIZE+CODE_SIZE][Q],cn[CODE_SIZE][A+2],vd[INFO_SIZE+CODE_SIZE],cd[CODE_SIZE];
  double trans[CODE_SIZE+INFO_SIZE],sig2;
  double llr[CODE_SIZE+INFO_SIZE],rate,ebn0,ebn0max,interval;
  int error,uplimit,num,j;
  
  srandom((unsigned)time(NULL));

#if NONSYS == 0
  //Rate calculation for systematic codes
  rate = (double)INFO_SIZE/((double)CODE_SIZE+(double)INFO_SIZE);
#else
  //Rate calculation for systematic codes
  rate = (double)INFO_SIZE/((double)CODE_SIZE);
#endif

  //show the resulting coding rate
  printf("coding rate: %lf\n",rate);
  
  ebn0     = atof(argv[1]);
  ebn0max  = atof(argv[2]);
  interval = atof(argv[3]);
  uplimit  = atoi(argv[4]);
  itemax   = atoi(argv[5]);

#if NONSYS == 1
  if (A != 1){
    printf("Fatal Error: Codes are catastrophic. A must be 1.");
    exit(1);
  }
#endif
  
  for (;;){ //EbN0 -> Max EbN0
    sig2 = 1.0/(rate * pow(10.0,ebn0/10.0)) * 0.5;
    error=0;
    
    for(num=0;num<uplimit;num++){
      info_gen(data); //information generation
      ra_random_interleaver(order); //interleaver for RA codes
      sys_ra_encoder(data,code,order); //encode      
      //transmission over AWGN channels
      for (i=0;i<CODE_SIZE+INFO_SIZE;i++){
	trans[i] = gauss(bpsk(code[i]),sqrt(sig2));
#if NONSYS == 1
	if (i > CODE_SIZE) trans[i] = 0.0; //puncturing systematic bits	
#endif
      }
      //setup tanner graph 
      generate_tanner_graph(order,vn,cn,vd,cd);

#if DETAIL == 1 //show details of generated codes
      show_graph_details(data,code,order,vd,cd,vn,cn);
#endif
      
      //calculate LLR at obs
      for (i=0;i<CODE_SIZE+INFO_SIZE;i++){
	llr[i] = llr_awgn(trans[i],sig2);
      }
      //sum-product decoding
      sum_product_decoder(llr,ans,vn,cn,vd,cd,itemax);

      //counting up # of errors
      for(i=0;i<INFO_SIZE;i++){
	if (data[i] != ans[i]) error++;
      }
    }
    if (error != 0){
      printf("%lf\t%e\n",ebn0,(double)error/((double)INFO_SIZE * (double)uplimit));
    }
    ebn0 += interval;
    if (ebn0 > ebn0max) break;
  }
}

void show_graph_details(int *data, int *code, int *order, int *vd, int *cd, int vn[][Q],int cn[][A+2]){
  int i,j;
  
  printf("data--\n");
  printintseq(data,INFO_SIZE);
  printf("code--\n");
  printintseq(code,INFO_SIZE+CODE_SIZE);
  printf("order---\n");
  printintseq(order,INFO_SIZE*Q);      
  printf("\novn---\n");
  for (i=0;i<CODE_SIZE;i++){
    for (j=0;j<vd[i];j++){
      printf("%d ",vn[i][j]);
    } printf("\n");
  }
  printf("\nhvn---\n");
  for (i=CODE_SIZE;i<INFO_SIZE+CODE_SIZE;i++){
    for (j=0;j<vd[i];j++){
      printf("%d ",vn[i][j]);
    } printf("\n");
  }
  printf("\ncn---\n");
  for (i=0;i<CODE_SIZE;i++){
    for (j=0;j<cd[i];j++){
      printf("%d ",cn[i][j]);
    } printf("\n");
  }
}

double llr_awgn(double y, double sig2){

  return(2.0*y/sig2);
  
}

void sum_product_decoder(double *llr,int *ans,int vn[][Q],int cn[][A+2],int *vd,int *cd, int itemax){
  double beta[CODE_SIZE+INFO_SIZE][Q]={};
  double alpha[CODE_SIZE][A+2]={};
  double alpha_sum,beta_prod,alpha_tmp[Q],beta_tmp[A+2];
  int i,k,j,tmp,n,focus,bitcheck,tmp_ans[CODE_SIZE+INFO_SIZE],l,flag=0;

  for (l=0;l<itemax;l++){
    if (l > 0){
      //check update---------------->
      for (i=0;i<CODE_SIZE;i++){
	for (n=0;n<cd[i];n++){ //target
	  beta_prod = 1.0;
	  for (k=0;k<cd[i];k++){
	    if (k!=n){
	      focus = cn[i][k]; //connected node
	      for (j=0;j<vd[focus];j++){
		if (vn[focus][j] == i){ tmp = j; break;}
	      }
	      beta_prod *= tanh(0.5*beta[focus][tmp]);
	    }
	  }
	  if (beta_prod > 0.9999999999999999){
	    alpha[i][n] = 2.0*atanh(0.9999999999999999);
	  }
	  else if (beta_prod < -0.9999999999999999){
	    alpha[i][n] = 2.0*atanh(-0.9999999999999999);
	  }else{
	    alpha[i][n] = 2.0*atanh(beta_prod);
	  }
	}
      }
    }
    
    //variable update ---------------->
    for (i=0;i<CODE_SIZE+INFO_SIZE;i++){
      alpha_sum=0.0;
      //calculating ith node
      for (k=0;k<vd[i];k++){ 
	//focus: connected check node
	focus = vn[i][k];
	for (j=0;j<cd[focus];j++){
	  if (cn[focus][j] == i){ tmp = j; break;}
	}
	alpha_sum += alpha[focus][tmp];
	alpha_tmp[k] = alpha[focus][tmp];
      }

      //decision - systematic case
      alpha_sum = llr[i] + alpha_sum;
      
      if (alpha_sum > 0.0){ tmp_ans[i] = 0;}
      else{ tmp_ans[i] = 1;}
      
      //update msg from ov to cn
      for (k=0;k<vd[i];k++){
	beta[i][k] = alpha_sum - alpha_tmp[k];	
      }
    }

    if (l>0){
      //check PC matrix ---------------->
      for (flag=1,i=0;i<CODE_SIZE;i++){
	bitcheck=0;
	for (j=0;j<cd[i];j++){
	  bitcheck = bitcheck ^ tmp_ans[cn[i][j]];
	}
	if (bitcheck != 0){ flag = 0; break;}
      }
      if (flag != 0) break;
    }
  }
  
  for (i=CODE_SIZE;i<CODE_SIZE+INFO_SIZE;i++){
    ans[i-CODE_SIZE] = tmp_ans[i];
  }

}

void generate_tanner_graph(int *order,int vn[][Q],int cn[][A+2],int *vd,int *cd){
  //maximum weight in column is Q
  //maximum weight in row is (A + 2)
  int i,k,tmp[INFO_SIZE*Q];

  //store degrees of v.n.
  for (i=0;i<CODE_SIZE+INFO_SIZE;i++){
    if (i < CODE_SIZE){
      if (i==CODE_SIZE-1){ vd[i] = 1; }
      else{
	vd[i] = 2;
      }
    }else{
      vd[i] = Q;
    }
  }
  //store degrees of c.n.
  for (i=0;i<CODE_SIZE;i++){
    if (i==0){ cd[i] = A+1; }
    else{
      cd[i] = A+2;
    }
  }
  
  //check nodes setup for ov
  for (i=0;i<CODE_SIZE;i++){
    for (k=0;k<cd[i]-A;k++){
	cn[i][k] = i - k;
    }
  }
  //check nodes setup for hv
  for (i=0;i<CODE_SIZE;i++){  
    for (k=cd[i]-A;k<cd[i];k++){
      cn[i][k] = (int)(order[A*i + (k-(cd[i]-A))]/Q) + CODE_SIZE;
    }
  }
  //ov nodes setup : (CODE_SIZE)
  for (i=0;i<CODE_SIZE;i++){
    for (k=0;k<vd[i];k++){
	vn[i][k] = i + k;
    }
  }
  //hv nodes setup : (INFO_SIZE)
  for (i=0;i<INFO_SIZE*Q;i++){
    vn[(int)(order[i]/Q) + CODE_SIZE][order[i]%Q] = (int)(i/A);
  }
  
}


double gauss(double av,double st)
{
  double u1,u2,u3,z;

     u1 = (double)random()/RAND_MAX;
     u3 = (double)sqrt(-2*(double)log(u1));
     u1 = (double)random()/RAND_MAX;
     z=(u3*(double)sin(2.0*M_PI*u1));
     return(av+st*z);
}

double bpsk(int b){

  //0 -> +1, 1 -> -1
  return(1.0 - 2.0*(double)b);
  
}

void printintseq(int *a,int size){
  int i;
  for (i=0;i<size;i++){
    printf("%d ",a[i]);   
  } printf("\n");
}

void printdoubleseq(double *a,int size){
  int i;
  for (i=0;i<size;i++){
    printf("%lf ",a[i]);
  } printf("\n");
}
		 
void info_gen(int *data){
  int i; 
  for(i=0;i<INFO_SIZE;i++){
    data[i] = (int)random() % 2;
  }
}

//systematic RA encoder (does not send original information bits).
void sys_ra_encoder(int *data, int *code, int *order){
  int i,j,k,tmp[Q*INFO_SIZE],comb,D;

  //repeat
  for (i=0;i<INFO_SIZE;i++){
    for (j=0;j<Q;j++){
      tmp[Q*i+j] = data[i];
    }
  }
  
  //interleave + combiner
  for (i=0;i<CODE_SIZE;i++){
    comb = 0;
    for (j=0;j<A;j++){
      comb = comb ^ tmp[order[i*A + j]];
    } code[i] = comb;
  }
  
  //accumulater
  D = 0; //initialize
  for (i=0;i<CODE_SIZE;i++){
    code[i] = D ^ code[i];
    D = code[i];
  }

  //add systematic part:
  for (i=CODE_SIZE;i<CODE_SIZE+INFO_SIZE;i++){
    code[i] = data[i-CODE_SIZE];
  }
  
}


//order will be a one-to-one mapping function 
void random_interleaver(int *order)
{
  int random1,dummy,i;
  
  for (i=0;i<Q*INFO_SIZE;i++){
    order[i] = i;
  }

  for (i=0;i<Q*INFO_SIZE;i++){
    for (;;){
      random1 = random() % INFO_SIZE*Q;
      if (random1 != i) break;
    }
    dummy = order[i];
    order[i] = order[random1];
    order[random1] = dummy;
  }
}

void ra_random_interleaver(int *order){

  struct member list;
  int i,j,tmp_intv[CODE_SIZE][A],shift,total_tmp,tmp,k;
  int tmp_order[INFO_SIZE*Q];
  
  list.next = NULL;
  for (k=0;k<A;k++){
    list_gen(CODE_SIZE,&list);
    total_tmp = CODE_SIZE;
    for (i=0;i<CODE_SIZE;i++){ //random interleave
      tmp = (int)random() % total_tmp;
      tmp_intv[i][k] = list_display(tmp,&list);
      list_del(tmp_intv[i][k],&list);
      total_tmp--;
    }
  }

  for (j=0,k=0;k<A;k++){
    for (i=0;i<CODE_SIZE;i++){
      tmp_order[j] = tmp_intv[i][k]*A + k;
      j++;
    }
  }

  for (i=0;i<INFO_SIZE*Q;i++){
    order[tmp_order[i]] = i;
  }
  
  list_release(&list);
  
}


void list_gen(int n, struct member *list){
  int i;
  for (i=0;i<n;i++){
    list_add(i,list);
  }
}

void list_add(int tmp, struct member *list){
  struct member *p;
  struct member *next;
  struct member *prev;

  //acquire a new memory region
  p = (struct member *)malloc(sizeof(struct member));
  if (p==NULL) {
  printf("out\n");
  }
  //plug given value into next pointer
  p->value = tmp;
  p->next  = NULL; //declear end point

  prev=list;
  
  for (next=list->next;next!=NULL;next=next->next){
    prev = next;
  }

  prev->next = p;
  
}

void list_del(int tmp, struct member *list){

  struct member *prev;
  struct member *p;

  prev = list;
  for (p=list->next;p!=NULL;p=p->next){
    if (p->value == tmp){
      if (p->next != NULL){
	prev->next=p->next;
	free(p);
	return;
      }
      prev->next = NULL;
      free(p);
      return;
    }
    prev = p;
  }

}


int list_display(int i, struct member *list){

  struct member *p;
  int j;

  if (list->next == NULL){
    printf("list is empty\n");
    exit(1);
  }

  for (j=0,p = list->next;p!=NULL&&j<i;j++){
    p=p->next;
  }

  return(p->value);

}

void list_all_display(struct member *list){

  struct member *p;

  if (list->next == NULL){
    printf("list is empty\n");
    exit(1);
  }

  for (p = list->next;p!=NULL;p=p->next){

    printf("%d ",p->value);

  } printf("\n");
}

void list_release(struct member *list){

  struct member *next;
  struct member *del;

  next = list->next;
  while (next){
    del=next;
    next=next->next;
    free(del);
  }
  list->next = NULL;

}
