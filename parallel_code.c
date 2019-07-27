#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#define CLK CLOCK_MONOTONIC

void integer_to_string(char ii[],int num)  //Converts integer to a string
{
	int a;
	int i=0;
	char ch;
	if(num==0)
	{
		ii[0]='0';
		ii[1]='\0';
		return;
	}
	while(num!=0)
	{
		a=num%10;
		num=num/10;
		ch=a+48;
//		printf("%c",ch);
		ii[i]=ch;
		i++;
	}
	ii[i]='\0';
	int len=i;
	char jj[len];
	for(i=0;i<len;i++)
	{
		jj[i]=ii[len-1-i];
	}
	for(i=0;i<len;i++)
	{
		ii[i]=jj[i];
	}
	ii[i]='\0';
//	printf("%s\n",ii);
}

struct timespec diff(struct timespec start, struct timespec end){  //Used to calculate the time taken to implement algorithm
    struct timespec temp;
    if((end.tv_nsec-start.tv_nsec)<0){
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    }
    else{
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}

typedef struct {
    unsigned char red,green,blue;
} PPMPixel;

typedef struct {
    int x, y;
    PPMPixel *data;
} PPMImage;

typedef struct {
    unsigned char gs;
} PPMPixelGS;


typedef struct {
    int x, y;
    PPMPixelGS *data;
} PPMImageGS;

//Two images are used: colour and greyscale images, images must be in .ppm format.

#define RGB_COMPONENT_COLOR 255


void writePPMGS(const char *filename, PPMImageGS *img) //Write a ppm image to a file
{
    FILE *fp;
    //open file for output
    fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Unable to open file '%s'\n", filename);
        exit(1);
    }

    //write the header file
    //image format
    fprintf(fp, "P5\n");



    //image size
    fprintf(fp, "%d %d\n",img->x,img->y);

    // rgb component depth
    fprintf(fp, "%d\n",RGB_COMPONENT_COLOR);

    // pixel data
    fwrite(img->data, img->x, img->y, fp);
    fclose(fp);
}


static PPMImage *readPPM(const char *filename) //read the ppm image from the file
{
    char buff[16];
    PPMImage *img;
    FILE *fp;
    int c, rgb_comp_color;
    fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Unable to open file '%s'\n", filename);
        exit(1);
    }

    //read image format
    if (!fgets(buff, sizeof(buff), fp)) {
        perror(filename);
        exit(1);
    }

    //check the image format
    if (buff[0] != 'P' || buff[1] != '6') {
        fprintf(stderr, "Invalid image format (must be 'P6')\n");
        exit(1);
    }

    //alloc memory form image
    img = (PPMImage *)malloc(sizeof(PPMImage));
    if (!img) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    //check for comments
    c = getc(fp);
    while (c == '#') {
        while (getc(fp) != '\n') ;
        c = getc(fp);
    }

    ungetc(c, fp);
    //read image size information
    if (fscanf(fp, "%d %d", &img->x, &img->y) != 2) {
        fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
        exit(1);
    }

    //read rgb component
    if (fscanf(fp, "%d", &rgb_comp_color) != 1) {
        fprintf(stderr, "Invalid rgb component (error loading '%s')\n", filename);
        exit(1);
    }

    //check rgb component depth
    if (rgb_comp_color!= RGB_COMPONENT_COLOR) {
        fprintf(stderr, "'%s' does not have 8-bits components\n", filename);
        exit(1);
    }

    while (fgetc(fp) != '\n') ;
    //memory allocation for pixel data
    img->data = (PPMPixel*)malloc(img->x * img->y * sizeof(PPMPixel));

    if (!img) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    //read pixel data from file
    if (fread(img->data, 3 * img->x, img->y, fp) != img->y) {
        fprintf(stderr, "Error loading image '%s'\n", filename);
        exit(1);
    }

    fclose(fp);
    return img;
}

void writePPM(const char *filename, PPMImage *img) //Write ppm image to file.
{
    FILE *fp;
    //open file for output
    fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Unable to open file '%s'\n", filename);
        exit(1);
    }

    //write the header file
    //image format
    fprintf(fp, "P6\n");

    //comments


    //image size
    fprintf(fp, "%d %d\n",img->x,img->y);

    // rgb component depth
    fprintf(fp, "%d\n",255);

    // pixel data
    fwrite(img->data, 3 * img->x, img->y, fp);
    fclose(fp);
}

void calculate_mean(int sum,int rows,int cols,double *means,double **a) //Calculate the mean image vector
{
	int i,j;
	for(i=0;i<cols;i++)
	{
		sum=0;
		for(j=0;j<rows;j++)
		{
			sum=sum+a[j][i];
		}
		means[i]=sum/rows;
	//	printf("%f\n",means[i]);
	}
}
void normalize_matrix(int rows,int cols, double **a,double **b,double *means)//Subtract the mean vector from each image vector
{
	int i,j;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			b[i][j]=a[i][j]-means[j];
		//	printf("%f ",b[i][j]);
		}
	//	printf("\n");
	//	printf("\n");
	}
}
void calculate_transpose(int rows,int cols,double **b,double **c) //Calculating the transpost of a matrix.
{
	int i,j;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			c[j][i]=b[i][j];
		}
	}
}
void initialize_eigenvectors(int cols,double **aaa)//Set all eigenvectors in the higher dimensional space to zero.
{
	int i,j;
	for(i=0;i<5;i++)
	{
		for(j=0;j<cols;j++)
		{
			aaa[i][j]=0;
		}
	}
}
void calculate_new_eigenvectors(int rows,int cols,double **aaa,double **eigenvectors,double **b)
{
	//Map the eigenvectors from a lower dimensional space to a higher dimensional space.
	int i,j,k;
	for(i=0;i<5;i++)
	{
		for(j=0;j<cols;j++)
		{
			for(k=0;k<rows;k++)
			{
				aaa[i][j]=aaa[i][j]+eigenvectors[i][k]*b[k][j];
			}
		//	printf("%f ",aaa[i][j]);
		}
	//	printf("\n");
	//	printf("\n");
	}
}
void initialize_resultant(int rows,double **resultant) //Intitialize the resultant where results of training phase are stored.
{
	int i,j;
	for(i=0;i<5;i++)
	{
		for(j=0;j<rows;j++)
		{
			resultant[i][j]=0;
		}
	}
}
void calculate_resultant(int rows,int cols,double **resultant,double **aaa,double **c)//Calculate results of training.
{
	int i,j,k;
	for(i=0;i<5;i++)
	{
		for(j=0;j<rows;j++)
		{
			for(k=0;k<cols;k++)
			{
				resultant[i][j]=resultant[i][j]+aaa[i][k]*c[k][j];
			}
		}
	}
}
void scale_resultant(int rows,double **resultant) //Training phase scaled by diving all values by 100000.
{
	int i,j;
	for(i=0;i<5;i++)
	{
		for(j=0;j<rows;j++)
		{
			resultant[i][j]=resultant[i][j]/100000;
	//		printf("%f ",resultant[i][j]);
		}
	//	printf("\n");
	}
//	printf("\n");
}
void initialize_new_resultant(double *resultant2) //Array to store test value initialized.
{
	int i;
	for(i=0;i<5;i++)
	{
		resultant2[i]=0;
	}
}
void calculate_new_resultant(int rrs,int cols,double **aaa,int aa[rrs][cols],double *resultant2) //Calculate resultant vector
{
	int i,j,k;
	for(i=0;i<5;i++)
	{
		for(j=0;j<1;j++)
		{
			for(k=0;k<cols;k++)
			{
				resultant2[i]=resultant2[i]+aaa[i][k]*aa[j][k];
			}
		}
	}
}
void find_suitable_matching_image(int no_of_images,double *resultant2,double **resultant,double *min_square_of_distances,int *image_no) //Compare with all images from the training set anf find which image
//the test image matches the best.
{
	int i,j;
	int image_number=0;
	double min=1000000000000000000;
	for(i=0;i<no_of_images;i++)
	{
		double distance=0;
		for(j=0;j<5;j++)
		{
			distance=distance+(int)((resultant2[j]-resultant[j][i])*(resultant2[j]-resultant[j][i]));
		//	printf("%f ",distance);
		}
		if(distance<min)
		{
			min=(int)distance;
			image_number=i;
		}
	}
	*min_square_of_distances=(int)min;
	*image_no=image_number;
//	printf("%d\n",*image_no);
}
int main(int argc, char* argv[]) {

    struct timespec start_e2e, end_e2e, start_alg, end_alg, e2e, alg;
    clock_gettime(CLK, &start_e2e);

    /* Check if enough command-line arguments are taken in. */
    /*******
	if(argc < 2){
        printf( "Usage: %s n p \n", argv[0] );
        return -1;
    }
*********/
    int N = atoi(argv[1]);                  /* size of input array */
//    int p = atoi(argv[2]);                  /* number of processors*/
    char *problem_name = "image_warping";
    char *approach_name = "data_division";

    FILE* outputFile;
    char outputFileName[100];
    sprintf(outputFileName, "output/%s_%s_%s_output.txt", problem_name, approach_name, argv[1]);

//    int number_of_threads = p;
//    omp_set_num_threads(p);
    char filename[1024];
    filename[0] ='\0';
    strcat(filename, argv[1]);
    strcat(filename, ".ppm");
    PPMImage *image;
//    image = (PPMImage *) malloc(sizeof(PPMImage));
	image = readPPM("image1.ppm");
//	free(image);
                    /* Start the algo timer */
   // PPMImageGS* x;// = change_image_warping(image);
    clock_gettime(CLK, &start_alg);
    //----------------------------------------Algorithm Here------------------------------------------
	
	double min_square_of_distances[1];
	double absolute_minimum=100000000;
	int image_no[1];
	int abs_image_no;
	int p,id;
	int ierr=MPI_Init(&argc,&argv);
	ierr=MPI_Comm_rank(MPI_COMM_WORLD,&id);
	ierr=MPI_Comm_size(MPI_COMM_WORLD,&p);	

		int i,j,idx;
		int total_images=N;
		int no_of_images=total_images/p;
		int no_of_rows=image->x;
		int no_of_columns=image->y;
		int cols=no_of_rows*no_of_columns;
		int rows=no_of_images;
		double **a=(double**)malloc(rows*sizeof(double*));
		for(i=0;i<rows;i++)
		{
			a[i]=(double*)malloc(cols*sizeof(double*));
		}
		double **b=(double**)malloc(rows*sizeof(double*));
		for(i=0;i<rows;i++)
		{
			b[i]=(double*)malloc(cols*sizeof(double*));
		}
		double *means=(double*)malloc(cols*sizeof(double));
		PPMImageGS *im2=(PPMImageGS *) malloc(sizeof(PPMImageGS));
		double sum=0;
		int no_of_images_read=0;
		for(i=id*no_of_images;i<no_of_images*(id+1);i++)
		{	
			char *file=(char*)malloc(1024*sizeof(char));
			file[0]='\0';
			strcat(file,"image");
			char *ii=(char*)malloc(5*sizeof(char));
			int num=0;
			integer_to_string(ii,i+1);
			strcat(file,ii);
			strcat(file,".ppm");
			image = readPPM(file);
			int no_of_rows=image->x;
			int no_of_columns=image->y;
			im2->x=no_of_rows;
			im2->y=no_of_columns;
			im2->data = (PPMPixelGS *) malloc(no_of_rows*no_of_columns*sizeof(PPMPixelGS));
			unsigned char img_values[no_of_rows*no_of_columns];
			for(j=0;j<cols;j++)
			{	
				im2->data[j].gs=(image->data[j].red+image->data[j].green+image->data[j].blue)/3;
				img_values[j]=im2->data[j].gs;
				a[no_of_images_read][j]=(double)img_values[j];
			//	a[no_of_images_read][j]=(double)i+j;
		//		printf("%d ",(int)a[no_of_images_read][j]);
			}
		//	printf("\n");
			no_of_images_read++;
		}
		calculate_mean(sum,rows,cols,means,a);
		normalize_matrix(rows,cols,a,b,means);
		double **c=(double**)malloc(cols*sizeof(double*));
		for(i=0;i<cols;i++)
		{
			c[i]=(double*)malloc(rows*sizeof(double*));
		}
		calculate_transpose(rows,cols,b,c);
		double **d=(double**)malloc(rows*sizeof(double*));
		for(i=0;i<rows;i++)
		{
			d[i]=(double*)malloc(rows*sizeof(double*));
		}
		double **eigenvectors=(double**)malloc(rows*sizeof(double*));
		for(i=0;i<rows;i++)
		{
			eigenvectors[i]=(double*)malloc(rows*sizeof(double*));
		}
		double *eigenvalues=(double*)malloc(rows*sizeof(double));
		int k;
		FILE *fp;
		char thread_id[5];
		integer_to_string(thread_id,id);
		char fn[30];
		fn[0]='\0';
		strcat(fn,"matrix");
		strcat(fn,thread_id);
		strcat(fn,".csv");
		fp=fopen(fn,"w+");
		for(i=0;i<rows;i++)
		{
			fprintf(fp,"Element,");
		}
		fprintf(fp,"\n");
		for(i=0;i<rows;i++)
		{
			for(j=0;j<rows;j++)
			{
				d[i][j]=0;
				for(k=0;k<cols;k++)
				{
					d[i][j]=d[i][j]+b[i][k]*c[k][j];
				}
				d[i][j]=d[i][j]/(cols-1);
				if(j!=cols-1)
				{
					fprintf(fp,"%f,",d[i][j]);
				}
				else
				{
					fprintf(fp,"%f",d[i][j]);
				}
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		fn[0]='\0';
		strcat(fn,"python eigs");
		strcat(fn,thread_id);
		strcat(fn,".py");
		system(fn);
	//	system(strcat("python eigs",strcat(thread_id,".py"));
	//	char *number;
		char number[15];
		size_t len;
		FILE *f1,*f2;
		fn[0]='\0';
		strcat(fn,"eigenvectors");
		strcat(fn,thread_id);
		strcat(fn,".txt");
		f1=fopen(fn,"r");
		fn[0]='\0';
		strcat(fn,"eigenvalues");
		strcat(fn,thread_id);
		strcat(fn,".txt");
		f2=fopen(fn,"r");
		for(i=0;i<rows;i++)
		{
			for(j=0;j<rows;j++)
			{
				int cc=0;
				while(1)
				{
					number[cc]=fgetc(f1);
					if(number[cc]==EOF || number[cc]=='\n')
					{
						number[cc]='\0';
						break;
					}
					else
					{
						cc++;
					}
				}
			//	getline(&number,&len,f1);
			//	printf("%s ",number);
				eigenvectors[i][j]=(double)atof(number);
		//		printf("%f ",eigenvectors[i][j]);
			}
		//	printf("\n");
		}
	//	printf("\n");
		for(i=0;i<rows;i++)
		{
			int cc=0;
			while(1)
			{
				number[cc]=fgetc(f1);
				if(number[cc]==EOF || number[cc]=='\n')
				{
					number[cc]='\0';
					break;
				}
				else
				{
					cc++;
				}
			}
	//		getline(&number,&len,f1);
			eigenvalues[i]=(double)atof(number);
	//		printf("%f ",eigenvalues[i]);
	//		printf("\n");
		}
		fclose(f1);
		fclose(f2);
	//	}
		for(i=0;i<rows;i++)
		{
			for(j=0;j<rows-1;j++)
			{
				if(eigenvalues[j]<eigenvalues[j+1])
				{
					double tmp=eigenvalues[j];
					eigenvalues[j]=eigenvalues[j+1];
					eigenvalues[j+1]=tmp;
					double tmp2[rows];
					for(k=0;k<rows;k++)
					{
						tmp2[k]=eigenvectors[j][k];
						eigenvectors[j][k]=eigenvectors[j+1][k];
						eigenvectors[j+1][k]=tmp2[k];
					}
				}
			}
		}
		double **aaa=(double**)malloc(5*sizeof(double*));
		for(i=0;i<5;i++)
		{
			aaa[i]=(double*)malloc(cols*sizeof(double*));
		}
		initialize_eigenvectors(cols,aaa);
		calculate_new_eigenvectors(rows,cols,aaa,eigenvectors,b);
		double **resultant=(double**)malloc(5*sizeof(double*));
		for(i=0;i<5;i++)
		{
			resultant[i]=(double*)malloc(rows*sizeof(double*));
		}
		initialize_resultant(rows,resultant);
		calculate_resultant(rows,cols,resultant,aaa,c);
		scale_resultant(rows,resultant);
		//Read the image to be tested./
		PPMImage *image2;
		char image_name[1024];
		image_name[0]='\0';
	//	strcat(image_name,argv[2]);
	//	strcat(image_name,"image7.ppm");
		image2 = readPPM("image7.ppm");
		int rrs=1;
		int aa[rrs][cols];
		unsigned char img_values2[no_of_rows*no_of_columns];
		PPMImageGS *im4=(PPMImageGS *) malloc(sizeof(PPMImageGS));
		im4->x=no_of_rows;
		im4->y=no_of_columns;
		im4->data = (PPMPixelGS *) malloc(no_of_rows*no_of_columns*sizeof(PPMPixelGS));
		
		for(i=0;i<cols;i++)
		{	
		//	idx=no_of_rows*i+j;
			im4->data[i].gs=(image2->data[i].red+image2->data[i].green+image2->data[i].blue)/3;
			img_values2[i]=im4->data[i].gs;
			aa[0][i]=(double)img_values2[i];
			aa[0][i]=aa[0][i]-means[i];
	//		printf("%d ",(int)aa[0][i]);
		}
		free(im4);
		double *resultant2=(double*)malloc(5*sizeof(double));
		//Compare with image and fine the suitable image, if no image matches, then print that no image was found.
		initialize_new_resultant(resultant2);
		calculate_new_resultant(rrs,cols,aaa,aa,resultant2);
		for(i=0;i<5;i++)
		{
			resultant2[i]=resultant2[i]/100000;
	//		printf("%f ",resultant2[i]);
		}
	//	printf("\n");
		find_suitable_matching_image(no_of_images,resultant2,resultant,&min_square_of_distances[0],&image_no[0]);
		/********
		if(id==0)
		{
			printf("%d ",image_no[0]);
		}
		*******/
		double msqd[p];
		int ino[p];
		if(id==0)
		{
			msqd[id]=min_square_of_distances[0];
			for(i=1;i<p;i++)
			{
				ierr=MPI_Recv(&msqd[i],1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			/*********
			for(i=0;i<p;i++)
			{
				printf("%f ",msqd[i]);
			}
		************/
		}		
		else
		{
			MPI_Send(min_square_of_distances,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		}
		if(id==0)
		{
			ino[id]=image_no[0];
			for(i=1;i<p;i++)
			{
				ierr=MPI_Recv(&ino[i],1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			/*******
			for(i=0;i<p;i++)
			{
				printf("%d ",ino[i]);
			}
			********/
		}
		else
		{
			MPI_Send(image_no,1,MPI_INT,0,0,MPI_COMM_WORLD);
		}
		if(id==0)
		{
			for(i=0;i<p;i++)
			{
				if(absolute_minimum>msqd[i])
				{
					absolute_minimum=msqd[i];
				abs_image_no=(rows*i)+ino[i];
				}
			}
		//	FILE *fppppp2;
		//	fppppp2=fopen("parallel_accuracy.csv","a");
		//	fprintf(fppppp2,"Value\n");
			if(absolute_minimum<1000)
			{
		//		if(abs_image_no+1==atoi(argv[2]))
		//		{
		//			fprintf(fppppp2,"%d\n",1);
		//		}
		//		else
		//		{
		//			fprintf(fppppp2,"%d\n",0);
		//		}
					printf("Image is of image %d.\n",abs_image_no);
			}
			else
			{
				printf("Image does not lie in the database.\n");
		//		fprintf(fppppp2,"%d\n",0);
			}
		//	fclose(fppppp2);
			clock_gettime(CLK, &end_alg);
			alg = diff(start_alg, end_alg);
			FILE *fpsssss;
			fpsssss=fopen("combined_logs.csv","a");
			fprintf(fpsssss,"%d,%d,%f\n",N,p,(double)alg.tv_sec+(double)alg.tv_nsec/1000000000);
			fclose(fpsssss);
		    printf("%s,%s,%d,%d,%ld,%ld\n", problem_name, approach_name, N, 0, alg.tv_sec, alg.tv_nsec);
		}
//	fclose(fppppp2);
	ierr=MPI_Finalize();
	
    //-----------------------------------------------------------------------------------------
     /* End the algo timer */
    char outputfilename[1024];
    outputfilename[0] ='\0';
    strcat(outputfilename, argv[1]);
    strcat(outputfilename, "_warped");
    strcat(outputfilename, ".ppm");
//    writePPMGS(outputfilename,im2);

    clock_gettime(CLK, &end_e2e);
    e2e = diff(start_e2e, end_e2e);
    
    
//    outputFile = fopen(outputFileName,"w");
	return 0;
}

