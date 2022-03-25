//SYLVAIN PERGAUD
//LAHCEN EL OUARDI
//TD1 TP1B

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define CAP_MAX 100
#define BUFSIZE 1024


//-----------------------------------------
//-----------------------------------------
//Rajout du .h car un seul fichier demander 
//-----------------------------------------
//-----------------------------------------

#ifndef __HULL__
#define __HULL__

//Creation de structure
struct vec {
	double x;
	double y;
};

struct vecset {
	struct vec *data;
	size_t size;
	size_t capacity;
};

//-------PARTIE 1------------------------------------------
//------PARTIE 1    -----------------
//----------------------------------------------------------------

/*
*Calcul le produit scalaire de deux vecteurs
*/
double dot(const struct vec *v1, const struct vec *v2);

/*
*Calcul le produit vectoriels de deux vecteurs definis pas 3 points
*/
double cross(const struct vec *p1, const struct vec *p2, const struct vec *p3);

/*
*Tournant a gauche en fonction du produit vectoriel 
*/
int is_left_turn(const struct vec *p1, const struct vec *p2, const struct vec *p3);

/*
*Permet de copier le vecteur v2 dans v1
*/
void copie(const struct vec *v1,  struct vec *v2);




//-------PARTIE 2------------------------------------------
//-------PARTIE 2   -----------------
//----------------------------------------------------------------



/*
*Create une structure vecset vide
*/
void vecset_create(struct vecset *self);

/*
*detruit un vecset
*/
void vecset_destroy(struct vecset *self);

/*
*ajoute un point à un ensemble de points
*/
void vecset_add(struct vecset *self, struct vec p);

/*
*fonction qui permet de comparer deux vecteur 
*/
int f_compare(const struct vec *p1, const struct vec *p2,char c);

/*
*fonction qui permet de comparer deux angle en fonction de leurs tangante  
*/
int f_compare_angle(const struct vec *p1, const struct vec *p2, const struct vec *p3);

/*
*fonction de comparaison
*/
typedef int (*comp_func_t)(const struct vec *p1, const struct vec *p2, const void *ctx);

/*
*Permet d'afficher un vecteur
*/
void vec_afficher(struct vec *v);

/*
*Permet d'afficher tout les vecteurs contenu dans le vecset et affiche aussi le size 
*/
void vecset_afficher(struct vecset *v);

/*
* renvoie le maximum d’un ensemble de points suivant une fonction de comparaison donnée 
*/
const struct vec *vecset_max(const struct vecset *self, comp_func_t func, const void *ctx);

/*
* renvoie le minimun d’un ensemble de points suivant une fonction de comparaison donnée 
*/
const struct vec *vecset_min(const struct vecset *self, comp_func_t func, const void *ctx);



/*void array_merge_sort_partial(struct vec *data,size_t i, size_t j, struct vec *tmp
				,comp_func_t func,const void *ctx);

void array_merge(struct vec *data, size_t i, size_t m, size_t j,struct vec *tmp,comp_func_t func,const void *ctx);

void vecset_sort(struct vecset *self, comp_func_t func,const void *ctx);
*/






void swap(const struct vecset *self, int i,int j);

/*
*fonction qui trie l’ensemble de points
suivant la fonction de comparaison donnée
*/
void sort(const struct vecset *self, comp_func_t func, const void *ctx);

/*
*fonction qui empile un élément
*/
void vecset_push(struct vecset *self, const struct vec p);

/*
*fonction qui dépile un élément
*/
void vecset_pop(struct vecset *self);

/*
*fonction qui renvoie le premier élément
de la pile
*/
const struct vec *vecset_top(const struct vecset *self);

/*
*fonction qui renvoie le second élément
de la pile
*/
const struct vec *vecset_second(const struct vecset *self);



int different(struct vec* p1, struct vec* p2);

//------------------------------------------
//-------PARTIE 3  JARVIS--------------------
//-----------------------------------------

/*
*fonction qui implémente la marche de JARVIS
*/
void jarvis_march(const struct vecset *in, struct vecset *out);


//-----------------------------------
//-------PARTIE 3 GRAHAM GRAHAM------
//-----------------------------------


/*
*fonction qui implémente le parcours de
Graham.
*/
void parcours_graham(struct vecset *in, struct vecset *out);

double norme(struct vec *v);

double distance(struct vec* X, struct vec* Y, struct vec* M);

/*
*fonction qui implémente le quickhull
*/

void quickhull(const struct vecset *in, struct vecset *out);

/*
*fonction qui implémente le findhull
*/
void findhull( struct vecset *in, struct vec *X,struct vec *Y, struct vecset *R);

#endif








//-------PARTIE IPLEMENTATION-------------
//--------PARTIE IPLEMENTATION-------------
//-----------------------------------------










double dot(const struct vec *v1, const struct vec *v2){
	return (v1->x*v2->x + v1->y*v2->y);
}

double cross(const struct vec *p1,const struct vec *p2, const struct vec *p3){
	return((p2->x-p1->x)*(p3->y-p1->y)-(p2->y-p1->y)*(p3->x-p1->x));
}


int is_left_turn(const struct vec *p1, const struct vec *p2, const struct vec *p3){
	
	return cross(p1,p2,p3)<0;
	

	
	

}
void copie(const struct vec *v1,  struct vec *v2){
	v2->x=v1->x;
	v2->y=v1->y;
	
}

/**** 

ENSEMBLE DE POINTS

*****/

void vecset_create(struct vecset *self){

	self->capacity = CAP_MAX;
	self->data=malloc(self->capacity*sizeof(struct vecset));
	self->size=0;


}

void vecset_destroy(struct vecset *self){
	free(self->data);
	self->data=NULL;

}

void vecset_add(struct vecset *self, struct vec p){
	if(self->capacity==self->size){
		self->capacity=self->capacity*2;
		struct vec *Mem=calloc(self->capacity,sizeof(struct vec));
		Mem=memcpy(Mem,self->data,self->size*sizeof(struct vec));
		free(self->data);
		self->data=Mem;
	}
	struct vec *p_tmp = malloc(sizeof(struct vec));
	p_tmp->x = p.x;
	p_tmp->y = p.y;
	self->data[self->size]=*p_tmp;
	self->size++;
	free(p_tmp);
	
}

int f_compare(const struct vec *p1, const struct vec *p2, char c){
	

	if(c=='a'){
		if(p1->x >p2->x){
			return -1;
		}else if(p1->x < p2->x){
			return 1;
		}
		else{
			if(p1->y > p2->y)
				return -1;
			else if(p1->y < p2->y)
				return 1;
			return 0;
		}
	}
	if(c=='o'){
		if(p1->y >p2->y){
			return -1;
		}else if(p1->y < p2->y){
			return 1;
		}
		else{
			if(p1->x > p2->x)
				return -1;
			else if(p1->x < p2->x)
				return 1;
			return 0;
		}
	}

}
int f_compare_angle(const struct vec *p1, const struct vec *p2, const struct vec *p3){

	double min = atan2(p2->y - p3->y, p2->x - p3->x);
	double tan = atan2(p1->y - p3->y, p1->x - p3->x);
	if(tan<min){
		return 1;
	} else if (tan>min){
		return -1;
	} else {
		return 0;
	}
}

void vec_afficher(struct vec *v){
	printf("%f %f\n",v->x,v->y); 
}

void vecset_afficher(struct vecset *v){
	printf("%zu\n",v->size);
	for(int i=0;i<v->size;i++){
		vec_afficher(&v->data[i]);
	
	}
	
}


const struct vec *vecset_max(const struct vecset *self,comp_func_t func, const void *ctx){
	if(self->size>0){
		struct vec *max=malloc(sizeof(struct vec));
		max->x=self->data[0].x;
		max->y=self->data[0].y;
		
		
		for(int i=1;i<self->size;i++){
			 if(func(&self->data[i],max,ctx)<0){
				max->x=self->data[i].x;
				max->y=self->data[i].y;
			 }
			
		}

		return max;
}return NULL;
}

const struct vec *vecset_min(const struct vecset *self, comp_func_t func, const void *ctx){


	if(self->size>0){
		struct vec *min=malloc(sizeof(struct vec));

		min->x=self->data[0].x;
		min->y=self->data[0].y;

	
		for(int i=1;i<self->size;i++){
			
			 if(func(&self->data[i],min,ctx)>0){
				min->x=self->data[i].x;
				min->y=self->data[i].y;
			 }
		
		}

		return min;
	}return NULL;

}






/*
void array_merge_sort_partial(struct vec *data,size_t i, size_t j, struct vec *tmp
				,comp_func_t func,const void *ctx) {
	if (j - i < 2) {
		return;
	}
	size_t m = (i + j) / 2;
	array_merge_sort_partial(data, i, m, tmp, func, ctx);
	array_merge_sort_partial(data, m, j, tmp, func, ctx);
	array_merge(data, i, m, j, tmp, func, ctx);
	memcpy(data + i, tmp + i, (j - i) * sizeof(struct vec));

}
void array_merge(struct vec *data, size_t i, size_t m, size_t j,struct vec *tmp,
			comp_func_t func,const void *ctx) {
	size_t a = i;
	size_t b = m;
	for (size_t k = i; k < j; ++k) {
		if (a < m && (b == j || func(&data[a], &data[b],ctx)<0)) {
			tmp[k].x = data[a].x;
			tmp[k].y = data[a].y;

			++a;
			} else {
				tmp[k].x = data[b].x;
				tmp[k].y = data[b].y;
				++b;
			}
	}
}



void vecset_sort(struct vecset *self, comp_func_t func,const void *ctx){
	if(self->size>0){
		struct vec *tmp = calloc(self->size, sizeof(struct vec));
		array_merge_sort_partial(self->data, 0, self->size, tmp,func,ctx);
		free(tmp);
	}
}
*/






void swap(const struct vecset *self, int i,int j){

	struct vec *tmp=malloc(sizeof(struct vec));
	copie(&self->data[i],tmp);
	copie(&self->data[j],&self->data[i]);
	copie(tmp,&self->data[j]);

}

void sort(const struct vecset *self, comp_func_t func, const void *ctx){
	for(int i=0;i<self->size;i++){
		int j=i;
		for(int k=j+1;k<self->size;k++){
			if(func(&self->data[k],&self->data[j],ctx)>0){
				j=k;
			}
		}

		swap(self,i,j);

		
		
	}


}
void vecset_push(struct vecset *self,const struct vec p){
	 vecset_add(self,p);
}

void vecset_pop(struct vecset *self){
	if(self->size>=1){
		self->size--;	
	}else{
		printf("Pile vide");
	}
	
}

const struct vec *vecset_top(const struct vecset *self){
	if(self->size>=1){
		return (&self->data[self->size-1]);
	}
	return NULL;
}


const struct vec *vecset_second(const struct vecset *self){
	if(self->size>=2){
		return(&self->data[self->size-2]);
		
	}
	return NULL;
}




//**********************************************************************
//************************************************************
//MARCHE DE JARVIS******************************************
//*****************************************************
//*****************************************






void jarvis_march(const struct vecset *in, struct vecset *out){

	comp_func_t f = f_compare;
	const struct vec *F = vecset_min(in,f,'a');
	struct vec C;
	copie(F,&C);
	int i=0;
	do {
	i++;
		vecset_push(out,C);
		struct vec N;
		copie(&in->data[i],&N);
		for(int i=0;i<in->size;i++){
			int is_left = is_left_turn(&C,&in->data[i],&N)>0;
			if(is_left){
				copie(&in->data[i],&N);
			}
		}

		copie(&N,&C);
	} while(F->x!=C.x && F->y!=C.y);
	
	free(F);
}



//******************
//ALGO de GRAHAM
//***********

void parcours_graham(struct vecset *in, struct vecset *out){

	comp_func_t f = f_compare;

	const struct vec *B = vecset_min(in,f,'o');

	comp_func_t f_a = f_compare_angle;
	sort(in,f_a,B);

	struct vec F = in->data[1];
	vecset_push(out,*B);
	vecset_push(out,F);
	struct vec *T=malloc(sizeof(struct vec));
	struct vec *S=malloc(sizeof(struct vec));
	for(int i=2;i<in->size;i++){	
		T=vecset_top(out);
		S=vecset_second(out);
		while(out->size >= 2 && is_left_turn(S,T,&in->data[i])){ 				
			vecset_pop(out);
			T=vecset_top(out);
			S=vecset_second(out);
		}
			
		vecset_push(out,in->data[i]);
	}
}
int different(struct vec *p1, struct vec* p2){
	return (p1->x != p2->x) || (p1->y != p2->y); 
}







//*****************************************
//************************************************
//ENVELOPPE RAPIDE******************************
//*****************************************
//*****************************************






double norme(struct vec *v){
	return sqrt(pow(v->x,2) + pow(v->y,2));

}

double distance(struct vec* X, struct vec* Y, struct vec* M){
	struct vec* vec1 = malloc(sizeof(struct vec));
	vec1->x = Y->x - X->x; vec1->y = Y->y - X->y;
	return fabs(cross(X,Y,M)/norme(vec1));
}  


void quickhull(const struct vecset *in, struct vecset *out){
	comp_func_t f = f_compare;
	struct vec *A = vecset_min(in,f,'a');



	struct vec *B = vecset_max(in,f,'a');


	
	struct vecset *S1=malloc(sizeof(struct vecset));
	vecset_create(S1);

	struct vecset *S2=malloc(sizeof(struct vecset));
	vecset_create(S2);

	for(int i=0;i<in->size;i++){
		if(different(&in->data[i],A) && different(&in->data[i],B)){
			if(is_left_turn(B,A,&in->data[i])){
				vecset_push(S1,in->data[i]);

			} else {
				vecset_push(S2,in->data[i]);

			}
		}
	}




	struct vecset* R1=malloc(sizeof(struct vecset));
	vecset_create(R1);
	//printf("S1:%zu\n",S1->size);
	findhull(S1,A,B,R1);
	struct vecset* R2=malloc(sizeof(struct vecset));
	vecset_create(R2);
	//printf("S2:%zu\n",S2->size);
	findhull(S2,B,A,R2);
	
	vecset_push(out,*A);
	for(int i=0;i<R1->size;i++){
		vecset_push(out,R1->data[i]);
	}
	vecset_push(out,*B);
	for(int i=0;i<R2->size;i++){
		vecset_push(out,R2->data[i]);
	}
}

void findhull( struct vecset *in, struct vec *X,struct vec *Y, struct vecset *R){
	if(in->size==0)
		return;

	double max = 0;
	struct vec M;
	copie(X,&M);

	for(int i=0;i<in->size;i++){
		if(different(X,&in->data[i]) && different(Y,&in->data[i])){

			double dist =distance(X,Y, &in->data[i]); 
				if(dist > max){
					copie(&in->data[i],&M);
					max = dist;
				}

			}
	}

	struct vecset* S1=malloc(sizeof(struct vecset));
	vecset_create(S1);
	struct vecset* S2=malloc(sizeof(struct vecset));
	vecset_create(S2);
	for(int i=0;i<in->size;i++){
		if(different(&in->data[i],&M)){
			if(is_left_turn(&M,X,&in->data[i])){
				vecset_push(S1,in->data[i]);
			}
			if(is_left_turn(Y,&M,&in->data[i])){
				vecset_push(S2,in->data[i]);
			}
		}
	}
	struct vecset* R1=malloc(sizeof(struct vecset));
	vecset_create(R1);
	findhull(S1,X,&M,R1);
	struct vecset* R2=malloc(sizeof(struct vecset));;
	vecset_create(R2);
	findhull(S2,&M,Y,R2);
	for(int i=0;i<R1->size;i++){
		vecset_push(R,R1->data[i]);
	}
	vecset_push(R,M);
	for(int i=0;i<R2->size;i++){
		vecset_push(R,R2->data[i]);
	}
}








//-----------------MAIN------------------------
//MAIN------------MAIN---------------------------
//---------------------MAIN--------------------
//-----------------------------------------
//-----------------------MAIN------------------


int main(){ 

	//Creation d'un vecset in
	struct vecset * in=malloc(sizeof(struct vecset));
	vecset_create(in);
	
	setbuf(stdout, NULL); 

	char buffer[BUFSIZE];
	fgets(buffer, BUFSIZE, stdin);
	size_t count = strtol(buffer, NULL, 10);
	for (size_t i = 0; i < count; ++i) {
	struct vec p;
	fgets(buffer, BUFSIZE, stdin);
	char *endptr = buffer;
	p.x = strtod(endptr, &endptr);
	p.y = strtod(endptr, &endptr);
	vecset_push(in,p);	//on rajoute les points dans le vecset in	
				
	}
	
	//Creation d'un vecset out
	struct vecset * out=malloc(sizeof(struct vecset));
	vecset_create(out);




/*Pour afficher JARVIS 
	-Décommenter le jarvis	
	-Compile avec " gcc *.c -lm -o hull"
	-execute avec "./hull-generator 100 | ./hull-viewer ./hull"
*/


	//jarvis_march(in,out);

		



/*Pour afficher GRAHAM 
	-Décommenter le graham
	-Compile avec " gcc *.c -lm -o hull"
	-execute avec "./hull-generator 100 | ./hull-viewer ./hull"
*/


	//parcours_graham(in, out);



/*Pour afficher l'enveloppe RAPIDE 
	-Décommenter quickhull
	-Compile avec " gcc *.c -lm -o hull"
	-execute avec "./hull-generator 100 | ./hull-viewer ./hull"

*/	
	
	quickhull(in,out);
	
	vecset_afficher(out);
	vecset_destroy(in);
	vecset_destroy(out);



	return 0;



}














