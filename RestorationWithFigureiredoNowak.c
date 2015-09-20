#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Function.h"


#define NAME_IMG_IN  "photograph"
#define NAME_IMG_OUT1 "photograph_original"
#define NAME_IMG_OUT2 "photograph_degraded" 
#define NAME_IMG_OUT3 "photograph_restaured" 
#define alpha 1 


int main(int argc,char** argv)
{
	int nb_iterations;
	int nbLevels=3;
  	int i,j,k,l,flag;
  	int length,width;
  	float var,sum1,sum2;
  	int size_filter;
	float** image;
	float** g;
	float** h;
	float** f;
	float** temp;
	float** temp1;
	float** haar;
	float** haar_inverse;
  		
  	printf("Input the size of low-pass filter : ");
  	scanf("%d",&size_filter);	

  	printf("Input the variance of noise : ");
  	scanf("%f",&var);
 
  	printf("Input the number of iterations (LANDWEBER): ");
  	scanf("%d",&nb_iterations);

	image=LoadImagePgm(NAME_IMG_IN, &length, &width);
	g=fmatrix_allocate_2d(length,width);
	f=fmatrix_allocate_2d(length,width);
	h=fmatrix_allocate_2d(length,width);
	temp=fmatrix_allocate_2d(length,width);
	temp1=fmatrix_allocate_2d(length,width);
	haar=fmatrix_allocate_2d(length,width);	
	haar_inverse=fmatrix_allocate_2d(length,width);

	for(i=0;i<length;i++){
		for(j=0;j<width;j++){
			g[i][j]=0.0;
			f[i][j]=0.0;
			h[i][j]=0.0;	
			temp[i][j]=0.0;
			temp1[i][j]=0.0;
			haar[i][j]=0.0;
			haar_inverse[i][j]=0.0;
		}
	}

  	/* add_gaussian_noise */
 	flou(h,size_filter,length,width);
	G(g,image,h,length,width);
	add_gaussian_noise(g,length,width,var);
	copy(f,g,length,width);

	k=1;
   	do{
		copy(temp1,f,length,width);
		sum1=0.0;
		sum2=0.0;
		for(i=0;i<length;i++){
			for(j=0;j<width;j++){
				sum1=sum1+SQUARE(image[i][j]-g[i][j]);
				sum2=sum2+SQUARE(image[i][j]-temp1[i][j]);
			}
		}
		printf("\n %d\t %f\t",k,10*log10(sum1/sum2));
	
		/* 1er etape : deconvolution with Landweber */
		for(l=0;l<nb_iterations;l++){
			update(temp1,temp,h,g,alpha,size_filter,length,width);
		}
		
		/* 2e etape : denoise with Haar coefficients */
		haar2D_complete(temp1,haar,nbLevels,length,width);
		for(i=0;i<length;i++){
			for(j=0;j<width;j++){
				if((i>length/pow(2,nbLevels)||j>width/pow(2,nbLevels))){
					if((SQUARE(haar[i][j])-3*SQUARE(var))>0){
						haar[i][j]=(SQUARE(haar[i][j])-3*SQUARE(var))/haar[i][j];
					}
					else{
						haar[i][j]=0.0;
					}
				}
				else{
					haar[i][j]=haar[i][j];	
				}
			}
		}
		ihaar2D_complete(haar,haar_inverse,nbLevels,length,width);
		flag=diff(haar_inverse,f,length,width);
		copy(f,haar_inverse,length,width);
		k++;
    	}while(k<=5);	
    	
	Recal(f,length,width);
		
    
	SaveImagePgm(NAME_IMG_OUT1,image,length,width);
	system("display photograph_original.pgm&");
	SaveImagePgm(NAME_IMG_OUT2,g,length,width);
	system("display photograph_degraded.pgm&");
	SaveImagePgm(NAME_IMG_OUT3,f,length,width);
	system("display photograph_restaured.pgm&");


	free_fmatrix_2d(image);  
	free_fmatrix_2d(f);
	free_fmatrix_2d(g);
	free_fmatrix_2d(h);
	free_fmatrix_2d(temp);
	free_fmatrix_2d(temp1);
	free_fmatrix_2d(haar);
	free_fmatrix_2d(haar_inverse);
  	
 	printf("\n Ending ... \n\n\n");
  	return 0; 	 
}
