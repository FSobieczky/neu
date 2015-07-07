/* /\******************************************************************************* */
/*  *                               neu.c */
/*  * */
/*  *  Purpose: pointing out 'untypical' features in digital image */
/*  *           Useful for narrowing down large amount of image data in pattern search */
/*  *  Method : Divide Image in Blocks N_x times N_y, look only at specific range of bits, then: */
/*  *           a.) count changes of greyvalues between neighboring pixels */
/*  *           b.) List the M most interesting such Blocks */
/*  *            */
/*  *  Input  : digital image in the ASCII-PGM Format */
/*  *  Output : List of M most `unusual' M blocks in form of little  */
/*                PGM-files, stored in directory called `PICTURES'  */
/*                and picture called file_orig.png,   where `file' is the  */
/*                input filename, showing highligthed patches indicating  */
/*                search results. */
/*              The program displays this file automatically (using the linux command `display'). */
/*              If the -v or -V option is invocted, then the bitplane picture is shown. */
/*              -V ONLY shows this bitplane picture, and does nothing else. */

/*  *  Usage  : neu <file> <N_x> <N_y> <begbit> <endbit> [number of results (default:30)] [-v|-V] */
/*  *  Help   : 'neu' (without arguments) explains usage */
/*  *           'neu -h' gives long help. */
/*  *  */
/*  *  Architecture :  main() includes (mainly): init_pixels(), init_regions(), main_loop() */
/*               main_loop() finally calls picturize_regions() (last part of code). */
/*  *\/ */

/* #include <stdlib.h> */
/* #include <stdio.h> */
/* #include <math.h> */
/* #include <string.h> */

/* #include "criterion.h" */

/* #define NO 1 */
/* #define YES 0 */

/* #define EXCLUSIVELY 2 */

/* #define MAX_B 11000  /\* largest picture width (broadness) *\/ */
/* #define MAX_H 11000  /\* largest picture height *\/ */

/* #define EPSILON 0.0000001 /\* Small number *\/ */

/* #define STANDARD_STRINGLENGTH 300     /\* used for filenames, system()-commands etc.*\/ */

/* #define REGIONNAMES "REGIONS"  /\* Directory name where the listed blocks are stored *\/ */
/* #define PICTUREDIR "PICTURES"  /\* Directory name where the listed blocks are stored *\/ */
/* #define PNGDIR "PNG"  /\* Directory name where the listed blocks are stored *\/ */

/* void init_regions(), main_loop(), throw_out(), generate_picture(); */
/* void make_pic_from_region(), picturize_regions(); */
/* void join(); */

/* float ratio_equal_pixels(); */
/* int n_equal_pixels(); */
/* double mean_of_region(); */
/* double var_of_region(); */
/* double mean_of_gradient(); */
/* double laplace_of_region(); */
/* double laplace_squared_of_region(); */
/* double sup_of_region(); */
/* int num_of_changes_of_region(); */
/* int doublecmpfunc(); */
/* int total_gradient_of_region(); */

/* static int pcompare(); */

/* /\* There are two 'layers': A pixel layer, and a region layer. First the picture *\/ */
/* /\* is loaded into the pixel layer (a linked list), and then a linked list of *\/ */
/* /\* regions is created which points to groups (=Blocks) of pixels). All the *\/ */
/* /\* calculations are done by going through the region list and computing the *\/ */
/* /\* respective values for each of the blocks. *\/ */

/* /\* The program is designed so that eventually the block structure can be changed into *\/ */
/* /\* arbitrary formed regions (not blocks). This feature is not yet implemented. *\/ */
/* /\* This is the reason for the linked-list technique used here (as opposed to simple arrays).*\/ */

/* char help_str[] = {"\n Usage: neu <file> <N_x> <N_y> <begbit> <endbit> [M] [-v|V|s]\n\n \ */
/* The `file' has to be in ASCII-PGM format (Magic Number P2) without comment line,\n\n \ */
/* `Nx' and `Ny' define the dimensions of the blocks into which the picture \n\ */
/*   is divided up\n\n\ */
/*  `begbit' and `endbit', both from {0,1,2,3,4,5,6,7} and begbit <= endbit, \n\ */
/*   is the range of the bits of the greyvalue image which will \n\ */
/*   actually be used,\n\n\ */
/*   `M' is the number of blocks from the first part of the list of 'interesting' blocks\n\ */
/*    that should be outputted. Default is 30.\n\ */
/* \n\n"};  */



/* struct pixelstruct           /\*  This is the struct for the Pixel-Level *\/ */
/* { */
/*   int c1; /\* The colors. Right now, only c1 is used = greyvalue *\/ */
/*   int c2; /\* These are all numbers between 0 and 255 *\/ */
/*   int c3; /\* ... *\/ */

/*   int px; /\* These are the integer X-Y coordinates of the pixel *\/ */
/*   int py; /\* in the image *\/ */

/*   struct pixelstruct *npixel; /\* next  *\/ */
/*   struct pixelstruct *fpixel; /\* former*\/ */
/*   struct pixelstruct *upixel; /\* upper *\/ */
/*   struct pixelstruct *lpixel; /\* lower --- NOT 'left' !!! *\/ */

/*   struct pixelstruct *exit_pixel; /\* this is set to NULL at boundary *\/ */

/*   struct regionstruct *parent;  /\* pointer the 'parenting' region *\/ */
/* } *pixel, *first_pixel; */

/* struct pixelstruct *goto_pixel();   */

/* struct regionstruct      /\* Definition of Pointer for Region-List *\/ */
/* { */
/*   struct pixelstruct *ipixel; */
/*   struct regionstruct *nregion;  */
/*   struct regionstruct *formregion; */
/*   int length; */
/* } *region, *firstregion, *lastregion; */


/*                   /\* Definition of Pointer for List of Neighbours *\/ */
/* struct neighbourstruct  */
/* {   */
/*   struct regionstruct *neireg;   */
/*   struct neighbourstruct *nextn; */
/*   struct neighbourstruct *formn; */
/* } *determine_neighbourhood(); */


/* int H, B;       /\* Global Height and Width (Width = "B'roadness"...) *\/ */
/* int Nx, Ny;                    /\* Broadness and Height of the Blocks *\/ */
/* int RN;                        /\* Remaining Pixel *\/ */

/* char *command_string, directory[STANDARD_STRINGLENGTH]; */
/* char pic_directory[STANDARD_STRINGLENGTH], png_directory[STANDARD_STRINGLENGTH]; */
/* char file[STANDARD_STRINGLENGTH], file_stem[STANDARD_STRINGLENGTH]; */
/* char orig_file[STANDARD_STRINGLENGTH], orig_file_temp[STANDARD_STRINGLENGTH]; */

/* int threshold=THRESHOLD3;  /\* Number of results to be displayed *\/ */
/* double overall_mean, overall_mean_gradient, overall_quotient, overall_var; */
/* double overall_laplace, overall_mean_laplace, overall_mean_sq_laplace; */
/* int overall_mean_num_changes=0; */

/* int verbosity=NO; */

/* main(int argc, char **argv) */
/* { */
/*   char command[STANDARD_STRINGLENGTH]; */
/*   char file[STANDARD_STRINGLENGTH];    */
/*                            /\* name of input file - 1st argument *\/ */

/*   int begbit, endbit;    /\* this interval of bits will be considered *\/ */
/*   int filelength;       /\* length of filename *\/ */

/*   if(argc==2) */
/*       if(strcmp(argv[1],"-h")==0) */
/* 	{ */
/* 	  printf("%s\n\n", help_str); */
/* 	  exit(1); */
/*         } */

/*   if(argc!=6 && argc!=7 && argc!=8) */
/*     { */
/*       printf("Usage: neu <file> <N_x> <N_y> <begbit> <endbit> [#(results)] [-v|V|s] \n\n       neu -h      (for help)\n\n");  */
/*       exit(-1); */
/*     } */

/*   strcpy(file, argv[1]);  /\* this is the input file name *\/ */
  
/*   Nx = atoi(argv[2]);     /\* horizontal Block dimension *\/  */
/*   Ny = atoi(argv[3]);     /\* vertical  Block dimension *\/ */

/*   begbit = atoi(argv[4]); /\* starting bit to be considered *\/ */
/*   endbit = atoi(argv[5]); /\* ending bit to be considered *\/ */

/*   if(argc==8) */
/*     { */
/*       if(strcmp(argv[7],"-v")==0) */
/* 	{ */
/* 	  verbosity = YES; */
/* 	  threshold=atof(argv[6]);   */
/* 	} */
/*       else if(strcmp(argv[7],"-V")==0) */
/* 	{ */
/* 	  verbosity = EXCLUSIVELY;	   */
/* 	  threshold=atof(argv[6]);   */
/* 	} */
/*       else if(strcmp(argv[7],"-s")==0) */
/* 	{ */
/* 	  verbosity = NO; */
/* 	  threshold=atof(argv[6]); */
/* 	} */
/*       else */
/* 	{ */
/* 	  printf("Usage: neu <file> <N_x> <N_y> <begbit> <endbit> [#(results)] [-v|V|s] \n\n");  */
/* 	  exit(-1); */
/* 	} */
/*     } */
/*   else if(argc==7) */
/*     { */
/*       if(strcmp(argv[6],"-v")==0) */
/* 	verbosity = YES; */
/*       else if(strcmp(argv[6],"-V")==0) */
/*         verbosity = EXCLUSIVELY; */
/*       else if(strcmp(argv[6],"-s")==0) */
/* 	verbosity = NO; */
/*       else */
/* 	threshold=atof(argv[6]); */
/*     } */
/*   else */
/*     { */
/*       threshold=THRESHOLD3; */
/*     } */

/*   filelength = (int)strlen(file); */
/*   if(file[filelength-4]=='.') */
/*     { */
/*       strncpy(file_stem, file, filelength-4); */
/*       if(verbosity!=NO) */
/* 	printf("-->%s\n", orig_file_temp); */
/*       sprintf(orig_file, "%s_orig.png", file_stem); */
/*       if(verbosity!=NO) */
/* 	printf("Decapitalized!\n"); */
/*     }   /\* makes sure the ending .PNG is in decaps *\/ */
/*   else */
/*     { */
/*       strcpy(file_stem, file); */
/*       sprintf(orig_file, "%s_orig.png", file_stem); */
/*     } /\* in this case, the file didn't have a .png/.PNG ending *\/ */

  
/*   /\* printf("Lengt of filename: %d; New Filename: %s\n", filelength, orig_file);  *\/ */

/*   /\* In init_pixels the picture is loaded into one single pixel-list *\/ */
/*   if(init_pixels(file, begbit, endbit)==NO)    */
/*     exit(-1);  /\* In this case the file argv[1] couldn't be opened. *\/ */

/*   init_regions();   /\* The parentregion in the pixels is set here. *\/ */
/*                     /\* and the a list of Nx times Ny sized *\/ */
/*   /\* pixel_lists is defined, to the beginning elements *\/ */
/*   /\* of which the regions of the regions list point *\/  */

/*   main_loop(begbit, endbit);  /\* Here, each region is checked *\/ */
/*                                     /\* and if found interesting sent *\/ */
/*                                     /\* to make_picture().  *\/ */
/*                                     /\* The results are stored in *\/ */
/*                                     /\* directory REGIONS *\/ */
/* } */

/* int init_pixels(char file[], int begbit, int endbit) */
/* { */
/*   int i,j, returnvalue, hoehe, breite, farben, colour, temp_color; */
/*   int max_greyvalue; */
/*   char version[STANDARD_STRINGLENGTH], rev[STANDARD_STRINGLENGTH]; */
/*   char magicn[STANDARD_STRINGLENGTH], temp[STANDARD_STRINGLENGTH]; */
/*   char qfile[STANDARD_STRINGLENGTH];  /\* name of quatized (output) file (~.qan) *\/ */
/*   char q_cmd[STANDARD_STRINGLENGTH+8];  */
/*   FILE *fp; */
/*   struct pixelstruct *pixel_new; */

/*   pixel = (struct pixelstruct*)malloc(sizeof(struct pixelstruct)+1); */

/*   returnvalue = YES;  /\* If something goes wrong, this will be `NO' *\/ */

/*   if((fp=fopen(file,"r"))==NULL) */
/*     { */
/*       printf("Couldn't open %s for reading !\n", file); */
/*       returnvalue = NO; */
/*     } */

/*   if( returnvalue == YES) */
/*     { */
/*       fscanf(fp,"%s",magicn); */

/*       fscanf(fp,"%d %d %d", &breite, &hoehe, &farben);    */

/*       B = breite;   /\* width of the picture *\/ */
/*       H = hoehe;    /\* height of the picture *\/ */

/*       if(verbosity!=NO) */
/* 	printf("The picture is %d x %d and has %d greyvalues.\n", B, H, farben); */

/*       for(i=0;i<H;i++) */
/* 	{ */
/* 	  for(j=0;j<B;j++) */
/* 	    { */
/* 	      fscanf(fp,"%d",&colour); */

/* 	      pixel_new=(struct pixelstruct*)malloc(sizeof(struct pixelstruct)+1); */

/* 	      /\* this somehow doesnt work:     *\/ */
/* 	      /\*     pixel_new->c1=turn_to_binary(colour,begbit,endbit); ... So:*\/ */
/* 	      temp_color = turn_to_binary(colour,begbit,endbit); */
/* 	      pixel_new->c1 = temp_color;    /\*  ... workaround.... *\/ */

/* 	      pixel_new->c3 = colour; /\* store the real thing here *\/ */
/* 	               /\* note: `color' is used here for greyvalue *\/ */

/* 	      pixel_new->px=j;   */
/* 	      pixel_new->py=i;  		 */
/* 	      pixel_new->npixel = (struct pixelstruct*)NULL; */
/* 	      pixel_new->exit_pixel = (struct pixelstruct*)NULL;	       */
/* 	      pixel_new->upixel = (struct pixelstruct*)NULL; */

/* 	      if(i!=0 || j!=0)  /\* if not the first *\/ */
/* 		{ */
/* 		  pixel_new->fpixel = pixel;  /\* former pixel *\/ */
/* 		  if(i==0)		     */
/* 		    pixel_new->lpixel = (struct pixelstruct*)NULL;  */
/* 		} */
/* 	      else   /\* if it IS the first: *\/ */
/* 		{ */
/* 		  pixel_new->fpixel = (struct pixelstruct*)NULL; */
/* 		  first_pixel = pixel_new; */
/* 		  if(verbosity!=NO) */
/* 		    printf("Here the  first pixel is stored:\n%d\n", pixel_new->c1); */
/* 		} */

/* 	      pixel->npixel=pixel_new;	 */
/* 	      pixel->exit_pixel=pixel_new; */
/* 	      pixel=pixel_new; */

/* 	      /\*	      printf("%d %d : %d  - ", i, j, pixel->c1); *\/ */
/* 	    } */
/* 	  if(verbosity!=NO) */
/* 	    printf("\n"); */
/* 	} */
    
/*       fclose(fp); */
/*      } */

/*   strcpy(qfile, file); */
/*   strcat(qfile, ".qan"); */

/*   if((fp=fopen(qfile,"w"))==NULL) */
/*     { */
/*       printf("Couldn't open %s for writing !\n", qfile); */
/*       returnvalue = NO; */
/*     } */
/*   else */
/*     { */
/*       pixel = first_pixel; */
/*       max_greyvalue = (int)pow(2, endbit); */
/*       fprintf(fp,"P2\n%d %d\n %d\n", B, H, max_greyvalue); */
/*       for(i=0;i<H;i++) */
/* 	{ */
/* 	  for(j=0;j<B;j++) */
/* 	    { */
/* 	      fprintf(fp, "%d ", pixel->c1); */
/* 	      pixel = pixel->npixel; */
/* 	    } */
/* 	  fprintf(fp, "\n"); */
/* 	} */
/*       fclose(fp); */
/*       if(verbosity!=NO) */
/* 	printf("Quantized picture written to file %s.\n", qfile); */
/*     } */

/*   if(verbosity==YES||verbosity==EXCLUSIVELY) */
/*     { */
/*       strcpy(q_cmd, "display "); */
/*       strcat(q_cmd, qfile); */
/*       strcat(q_cmd, " &"); */
/*       system(q_cmd);  */
/*     } */
    
/*   if(verbosity==EXCLUSIVELY) /\* only show the .qan file with the chosen bitplanes*\/ */
/*     exit(0); */

/*   if(verbosity!=NO){ */
/*     printf("\n----------------------------------------\n"); */
/*     printf("Region '%s': Height times Width = %d x %d\n",file,H,B); */
/*   } */
    
/*   return returnvalue; */
/* } */

/* void init_regions() */
/* { */
/*   int i,j,k,l; */
/*   struct pixelstruct *cpixel, *ptemp, *terminal_pixel; */
/*   struct regionstruct *current, *lastcurrent, *nextcurrent;   */

/*   terminal_pixel=(struct pixelstruct*)malloc(sizeof(struct pixelstruct)+1); */
/*   terminal_pixel->c1=-1; */

/*   current = (struct regionstruct*)NULL; */
/*   if(verbosity!=NO) */
/*     printf("Now the regions:\n"); */

/*   RN = 0; */

/*   for(i=0; i<H-Ny; i=i+Ny) */
/*     for(j=0; j<B-Nx; j=j+Nx) */
/*       {         */
/* 	nextcurrent = (struct regionstruct*)malloc(sizeof(struct regionstruct)+1); */
/* 	RN++; */

/* 	if(i==0 && j==0) */
/* 	  firstregion =nextcurrent; */

/* 	cpixel = goto_pixel(j, i, B, H, first_pixel); */

/* 	nextcurrent->ipixel = cpixel; */

/* 	for(k=0;k<Ny;k++) */
/* 	  for(l=0;l<Nx;l++) */
/* 	    {	       */
/* 	      if(k!=0 || l!=0)  /\* if not the first pixel in that region *\/ */
/* 		{		   */
/* 		  cpixel->exit_pixel = goto_pixel(j+l, i+k, B, H, first_pixel); */
/* 		  cpixel = cpixel->exit_pixel;		  */

/* 	          cpixel->parent = nextcurrent; */
/* 		} */
/* 	    }	 */

/* 	cpixel->exit_pixel=(struct pixelstruct*)NULL; */
/* 	if(current!=(struct regionstruct*)NULL) */
/* 	    current->nregion = nextcurrent; */
/* 	nextcurrent->formregion=current; */
/* 	nextcurrent->nregion=(struct regionstruct*)NULL; */
/* 	lastcurrent = current; */
/* 	current = nextcurrent; */
	
/* 	if(verbosity!=NO) */
/* 	  printf("(%d %d) in (%d %d)\n", j, i, B, H); */
/*       } */

/*   lastregion=current; */

/*   firstregion->formregion =  (struct regionstruct*)NULL; */

/*   if(verbosity!=NO) */
/*     printf("\nInitialization of regions done:   We have %d   %d x %d  'Blocks'.\n",RN, Nx, Ny); */

/* } */

/* struct pixelstruct *goto_pixel(x, y, B, H, initial_pixel) */
/*      int x, y, B, H; struct pixelstruct *initial_pixel; */
/* { */
/*   int i, j; */
/*   struct pixelstruct *c_pixel; */

/*   c_pixel = initial_pixel; */

/*   for(i=1;i<y;i++) */
    for(j=1;j<=B;j++)	
      c_pixel = c_pixel->npixel;

  for(j=1;j<=x;j++)
    c_pixel = c_pixel->npixel;

  return c_pixel;
}

void main_loop(int begbit, int endbit)
{
  char statfile[STANDARD_STRINGLENGTH];  
      /\* name of stat-results (output)file(~.stat): use with R *\/  

  int i,j,k,sum, inc=0, count=0, count2=0, l=0, N_Blocks, fac;
  int freq_liste[(int)H*B/(Nx*Ny)+1], n_inc;
  int ok_to_picturize[(int)H*B/(Nx*Ny)+1];
  int region_num_changes[(int)H*B/(Nx*Ny)+1];
  int totalgrad[(int)H*B/(Nx*Ny)+1];  
  double p[(int)H*B/(Nx*Ny)+1], Entropie, total; 
   double region_mean[(int)H*B/(Nx*Ny)+1], region_gradient[(int)H*B/(Nx*Ny)+1]; */
  double region_var[(int)H*B/(Nx*Ny)+1], region_quotient[(int)H*B/(Nx*Ny)+1];
  double region_laplace[(int)H*B/(Nx*Ny)+1], region_sup[(int)H*B/(Nx*Ny)+1];
  double totgrad[(int)H*B/(Nx*Ny)+1];
  double region_sq_laplace[(int)H*B/(Nx*Ny)+1], Performance[(int)H*B/(Nx*Ny)+1];
  double Copy_of_Performance[(int)H*B/(Nx*Ny)+1];
  /* the `overall_...' variables for the stats are global (like region_mean[]...) */
  double threshold2=THRESHOLD2;   /* this parameter doesn't change so much */
  struct regionstruct *current, *temp, *temp2, *list[(int)H*B/(Nx*Ny)+1];;
  struct neighbourstruct *currentneigh, tempneigh;
  float f_sum = 0.0, f_inc, count3;
  FILE *fp;

  if(verbosity!=NO){
    if(begbit==endbit)
      printf("Considered: Bit %d.\n",endbit);
    else
      printf("Considered: Bits %d bis %d.\n", begbit, endbit);
  }

  N_Blocks = (int)H*B/(Nx*Ny);

  current = firstregion;     /* begin of the linked list of regions */

  currentneigh = (struct neighbourstruct*)malloc(sizeof(struct neighbourstruct)+1);
  currentneigh->neireg = firstregion;
  currentneigh->formn = NULL;

  while(current->nregion!=NULL)   /* Loops through regions with `current' */
    {
      if(verbosity!=NO)
	printf("Inside Region #%d with initial pixel (%d, %d)\n", count, current->ipixel->px, current->ipixel->py);
      temp = firstregion->nregion;
      sum = 0;
      f_sum = 0.0;

      region_mean[count] = mean_of_region(current);
      region_var[count] = var_of_region(current, region_mean[count]);
      region_gradient[count] = mean_of_gradient(current);
      region_laplace[count] = laplace_of_region(current);
      region_sq_laplace[count] = laplace_squared_of_region(current);
      totalgrad[count] = total_gradient_of_region(current);
      region_sup[count] = sup_of_region(current);
      region_quotient[count] = region_sq_laplace[count]/region_var[count];
      region_num_changes[count] = num_of_changes_of_region(current);

      list[count] = current;

    /*  ... this is where the part goes from archive.c that concerns the entropy .... not implemented here */

      if(verbosity!=NO)
	printf("n=%d: rms[n]=%5.3f,  region_gradient[n]=%5.3f,  quotient[n]=%5.3f\n", count, sqrt(region_var[count]), region_gradient[count], region_quotient[count]);

      current = current->nregion;
      count++;
   }

  strcpy(statfile, file);
  strcat(statfile, ".stat");

  /* Statistics of Regions' means, mean gradients, and variances - and their quotients*/
  
  overall_mean = 0.0;
  overall_var = 0.0;
  overall_laplace = 0.0;
  overall_mean_gradient = 0.0;
  overall_mean_laplace = 0.0;
  overall_mean_sq_laplace = 0.0;
  overall_quotient = 0.0;
  overall_mean_num_changes = 0;
  for(i=0;i<count;i++)
    {
      overall_mean += region_mean[i];
      overall_var += region_var[i];
      overall_mean_gradient += region_gradient[i];
      overall_quotient += region_quotient[i];
      overall_mean_laplace += region_laplace[i];
      overall_mean_sq_laplace += region_sq_laplace[i];
      overall_mean_num_changes += region_num_changes[i];
    }                   /* calculate mean over all blocks of these quantities */
  overall_mean = overall_mean/count;
  overall_var = (overall_var-overall_mean*overall_mean)/count;
  overall_mean_gradient = overall_mean_gradient/count;
  overall_quotient = overall_quotient/count;
  overall_mean_laplace = overall_mean_laplace/count;
  overall_mean_sq_laplace= overall_mean_sq_laplace/count;
  overall_mean_num_changes= overall_mean_num_changes/count;

  if(verbosity!=NO)
    printf("Mean: %9.5f - RMS: %9.5f -  Mean-Gradient: %9.5f - Mean quotient: %9.5f - Overall quotient: %9.5f\n", overall_mean, sqrt(overall_var), overall_mean_gradient, overall_quotient, overall_mean_gradient/sqrt(overall_var));

  /* getchar(); */

  if((fp=fopen(statfile,"w"))==NULL)
    {
      printf("Couldn't open %s for writing !\n", statfile);
    }
  else
    {
      fprintf(fp,"# Stats-Data for a single run of 'neu'.\n");
      fprintf(fp,"# 1.n, 2.Overall Mean, 3.Overall Var, 4.Mean, 5.Variance, 6.Gradient, 7.Laplace, 8.Lap_sq, 9.Num_Ch., 10.Quotient (8./5.)\n");
      for(i=1;+i<=count;i++)
	{
	  fprintf(fp, "%d, %9.6f, %9.6f, %9.6f, %9.6f, %9.6f, %9.6f, %9.6f, %d %9.6f\n ", i, overall_mean, overall_var, region_mean[i], region_var[i],region_gradient[i], region_laplace[i], totalgrad[i], region_num_changes[i], region_quotient[i]);
	}
      fclose(fp);
    }

  for(i=1;i<=count;i++)
    {
      Performance[i] = 1.0*G(i,overall_mean,overall_var,region_mean[i],region_var[i],region_gradient[i],region_laplace[i],totalgrad[i],region_num_changes[i],threshold,threshold2);
      Copy_of_Performance[i] = Performance[i];
    }

  qsort(Performance, count, sizeof(double), doublecmpfunc);
    /* Performance is sorted, Copy_of_Performance not */
    /* in this way, it is possible to find out for which i the */
    /* i'th region has high performance */

  /* for(i=1;i<=count;i++)
    {
      printf("Performance[%d]: %9.6f - Original[%d]: %9.6f\n", i,Performance[i],i,Copy_of_Performance[i]);
    }
  */

  for(i=1;i<=count;i++)
    {
      if(Copy_of_Performance[i]>=Performance[count-threshold])
	{
	  ok_to_picturize[i]=YES;
	  if(verbosity!=NO)
	    printf("Performance[%d]: %9.6f - Original[%d]: %9.6f\n", j,Performance[j],i,Copy_of_Performance[i]);
	  }
      else
	{
	  ok_to_picturize[i]=NO;
	}
    }

  picturize_regions(firstregion, begbit, endbit, ok_to_picturize, count, list);
  if(verbosity!=NO)
    system("./show");
}

/* comp_bits() vergleicht current and temp und sendet 1 bei Gleichheit, sonst, 0 */

int comp_bits(current,temp)
struct regionstruct *current, *temp;
{
  int value=1, go_on=1;
  int k, p1, p2, colour1, colour2;
  struct pixelstruct *pixel1, *pixel2;

  pixel1 = current->ipixel;
  pixel2 = temp->ipixel;

  while(go_on==1)
    {
      colour1 = pixel1->c1;
      colour2 = pixel2->c1;

      if(colour1!=colour2)
	value = 0;

      if(pixel1->exit_pixel==NULL || pixel1->exit_pixel==NULL)
	go_on=0;
      else
	{
	  pixel1 = pixel1->exit_pixel;
	  pixel2 = pixel2->exit_pixel;
	}
	
    }
  if(verbosity!=NO)
    printf("Everything ok!\n");

  /*  printf("Value returned: %d\n",value); */

  return value;
}

/* ratio_equal_pixels() vergleicht current and temp und sendet 1 bei
   Gleichheit, sonst, 0 */

/* ratio_equal_pixels() compares current and temp and sends 1 if they are equal, otherwise it sends 0. */

float ratio_equal_pixels(current,temp) 
struct regionstruct *current, *temp;
{
  float value;
  int n_equals=0, go_on=1;
  int k, p1, p2, colour1, colour2, n_pixels=0;
  struct pixelstruct *pixel1, *pixel2;

  pixel1 = current->ipixel;
  pixel2 = temp->ipixel;

  while(go_on==1)
    {      
      colour1 = pixel1->c1;
      colour2 = pixel2->c1;

      n_pixels++;
      if(colour1==colour2)
	n_equals += 1;

      if(pixel1->exit_pixel==NULL || pixel1->exit_pixel==NULL)
	go_on=0;
      else
	{
	  pixel1 = pixel1->exit_pixel;
	  pixel2 = pixel2->exit_pixel;
	}
    }

  value = (float)n_equals/n_pixels;

  if(verbosity!=NO)
    printf("Everything ok! (Ratio = %5.3f\n", value);

  return value;
}

/* n_equal_pixels() vergleicht current and temp und sendet Anzahl gleicher Pixel */

int n_equal_pixels(current,temp) 
struct regionstruct *current, *temp;
{
  int value;
  int n_equals=0, go_on=1;
  int k, p1, p2, colour1, colour2, n_pixels=0;
  struct pixelstruct *pixel1, *pixel2;

  pixel1 = current->ipixel;
  pixel2 = temp->ipixel;

  while(go_on==1)
    {      
      colour1 = pixel1->c1;
      colour2 = pixel2->c1;

      n_pixels++;
      if(colour1==colour2)
	n_equals += 1;

      if(pixel1->exit_pixel==NULL || pixel1->exit_pixel==NULL)
	go_on=0;
      else
	{
	  pixel1 = pixel1->exit_pixel;
	  pixel2 = pixel2->exit_pixel;
	}
    }

  value = n_equals;

  if(verbosity!=NO)
    printf("Everything ok! (# of equal pixels = %d)\n", value);

  return value;
}


/* wirft ein Listenelement raus -  throws out one list-element */

void throw_out(struct regionstruct *temp)
{
  if(temp->nregion!=NULL)
    {
      temp->formregion->nregion = temp->nregion;
      temp->nregion->formregion = temp->formregion;
    }
  else
    temp->formregion->nregion = NULL;
      
  free(temp);  

  RN--;
}

/* here the Bits are removed that are not */
/*   in the interval {begbit, ...., endbit} */

turn_to_binary(int d, int begbit, int endbit)
{
  int m, value = 0;
  
  for(m=7;m>=0;m--)
    {
      if((int)pow(2,m)<=d)
	{
	  if(m>=begbit && m<=endbit)
	    value += pow(2,m);
	  d = d - (int)pow(2,m);
	}
    }

  return value;
}


       /* Greyvalue-Mean of the pixels in the region: */
double mean_of_region(struct regionstruct *current)
{
  struct pixelstruct *pixel;
  int greyvalues=0, n_pixels=0;
  double value;

  pixel = current->ipixel;

  while(pixel!=NULL)
    {
      greyvalues += pixel->c1;
      pixel = pixel->exit_pixel;
      n_pixels++;
    }

  value = (double)greyvalues/n_pixels;

  return value;
}

/* ... and corresponding variance: */
double var_of_region(struct regionstruct *current, double mean)
{
  struct pixelstruct *pixel;
  int greyvalues=0, n_pixels=0, temp_value;
  double value;

  pixel = current->ipixel;

  while(pixel!=NULL)
    {
      temp_value = pixel->c1;
      greyvalues += temp_value*temp_value;
      pixel = pixel->exit_pixel;
      n_pixels++;
    }

  value = (double)(greyvalues/n_pixels - mean*mean);

  return value;
}

      /* this spits out the mean of the greyvalue-differences */
      /* of consecutive pixels in the list (corresponds to gradient */
      /* in X-direction). */
double mean_of_gradient(struct regionstruct *current)
{
  struct pixelstruct *pixel, *ipixel;
  int gradient_value=0, n_pixels=0, go_on=1;
  double value;

  ipixel = current->ipixel;
  pixel = ipixel->exit_pixel;

  while(pixel!=NULL)
    {
      gradient_value += abs(pixel->c1 - (pixel->fpixel)->c1);      
      pixel = pixel->exit_pixel;  
      n_pixels++;
     }

  value = (double)gradient_value/n_pixels;
  
  return value;
}

      /* this produces the mean of the consecutive differences of the differences */
      /* of consecutive pixels in the pixel-list of the given region */
double laplace_of_region(struct regionstruct *current)
{
  struct pixelstruct *pixel, *ipixel;
  int laplace_value=0, n_pixels=0, go_on=1;
  double value;

  ipixel = current->ipixel;
  pixel = ipixel->exit_pixel; /* recall: exit_pixel is n_pixel, or NULL if at */
                              /*         end of region */

  while(pixel->exit_pixel!=NULL)
    {
      laplace_value += abs((pixel->exit_pixel->c1 - pixel->c1)-(pixel->c1 - pixel->fpixel->c1));
      pixel = pixel->exit_pixel;  
      n_pixels++;
     }

  value = (double)laplace_value/n_pixels;
  
  return value;
}

      /* this produces the mean of the 2nd power of consecutive differences  */
      /* of the differences of consecutive pixels in the list */
double laplace_squared_of_region(struct regionstruct *current)
{
  struct pixelstruct *pixel, *ipixel;
  int laplace_value=0, n_pixels=0, go_on=1;
  double value;

  ipixel = current->ipixel;
  pixel = ipixel->exit_pixel; /* recall: exit_pixel is n_pixel, or NULL if at */
                              /*         end of region */

  /*  while(pixel->exit_pixel!=NULL)*/
  while(pixel!=NULL)
    {
      laplace_value += pow((pixel->npixel->c1 - pixel->c1)-(pixel->c1 - pixel->fpixel->c1), 2);
      pixel = pixel->exit_pixel;  
      n_pixels++;
     }

  value = (double)sqrt(laplace_value)/n_pixels;
  
  return value;
}


      /* this gives the difference of the maximum and minimum of the */
      /* c1-values (greyvalues) occuring in the region */
double sup_of_region(struct regionstruct *current)
{
  struct pixelstruct *pixel, *ipixel;
  int laplace_value=0, n_pixels=0, go_on=1;
  int max_value=0, min_value=255;
  double value;

  ipixel = current->ipixel;
  pixel = ipixel; /* recall: exit_pixel is n_pixel, or NULL if at */
                              /*         end of region */

  while(pixel!=NULL)
    {
      if(pixel->c1 > max_value)
	max_value = pixel->c1;
      if(pixel->c1 < min_value)
	min_value = pixel->c1;

      pixel = pixel->exit_pixel;  
     }

  value = (double)(max_value-min_value)/255.0;
  
  return value;
}

int num_of_changes_of_region(struct regionstruct *current)
{
  struct pixelstruct *pixel, *ipixel;
  int num_changes=0, n_pixels=0, go_on=1;
  int max_value=0, min_value=255;
  double value;

  ipixel = current->ipixel;
  pixel = ipixel; /* recall: exit_pixel is n_pixel, or NULL if at
                              /*         end of region */

  while(pixel->exit_pixel!=NULL)
    {
      if(pixel->c1 != pixel->npixel->c1)
	num_changes++;

      pixel = pixel->exit_pixel;  
    }  

  return num_changes;
}

int total_gradient_of_region(struct regionstruct *current)
{
  struct pixelstruct *pixel, *ipixel;
  int tgrad=0, n_pixels=0, go_on=1;
  int max_value=0, min_value=255;
  int gradtemp;
  double value;

  ipixel = current->ipixel;
  pixel = ipixel; /* recall: exit_pixel is n_pixel, or NULL if at */
                              /*         end of region */

  while(pixel->exit_pixel!=NULL)
    {
      tgrad+= (pixel->c1 - pixel->npixel->c1);
      pixel = pixel->exit_pixel;  
    }  

  return tgrad;
}


int doublecmpfunc (const void * a, const void * b)
{
   return ( *(double*)a - *(double*)b );
}


/*-------               Picturize - Section            --------*/


void picturize_regions(struct regionstruct *first_cur, int begbit, int endbit, int ok_to_picturize[], int length_of_region_list, struct regionstruct *list[])
{
  int i,j,colour, x, y, n, go_on=1, len_region, count=0;
  struct regionstruct *current;
  struct pixelstruct *cpix;
  char reg_fname[STANDARD_STRINGLENGTH];     
  char squares_fname[STANDARD_STRINGLENGTH]; /* this has the little boxes painted in it*/
  char ofile_cmd[STANDARD_STRINGLENGTH];     
  FILE *fp, *fp1;
  int Original_Picture[B*H];
  int highlighted_color;

  command_string = (char*)malloc(sizeof(char)*RN*50+1); 
                  /* command_string is a global variable */

  current = first_cur;  /* beginning of list */

  for(i=0;i<length_of_region_list;i++)
    {
      cpix=list[i]->ipixel;
      while(cpix!=NULL)
	{
	  x=(int)cpix->px;
	  y=(int)cpix->py;
	  if(ok_to_picturize[i]==NO)
	    Original_Picture[y*H+x]=(int)cpix->c3;
	  else
	    {
	      highlighted_color=cpix->c3 + 30;
	      if(highlighted_color>255)
		highlighted_color=255;
	      Original_Picture[y*H+x]=highlighted_color;
	    }
	  cpix=cpix->exit_pixel;
	}      
    } 

  if((fp1=fopen(orig_file,"w"))==NULL)
    {
      printf("Couldn't open %s for writing !\n", orig_file);
      printf("No picture with boxes inside will be produced!\n");
    }
  else
    {
      fprintf(fp1,"P2\n");
      fprintf(fp1,"%d %d \n%d\n", B,H, 255);

      for(i=0;i<H;i++)
	{
	  for(j=0;j<B;j++)
	    fprintf(fp1,"%d ", Original_Picture[j+H*i]);
	  fprintf(fp1, "\n");
	}
      fclose(fp1);
    }

  if(verbosity!=NO){
    sprintf(ofile_cmd,"display %s &", orig_file);
    system(ofile_cmd);
  }

  current = first_cur;  /* beginning of list */
  n = 0;  /* Counts number of pictures made */

  while(current!=NULL)
    {
      len_region = length_region(current);
      /* printf("Length of the region: %d\n", len_region);*/
      if(ok_to_picturize[count]==YES) 
	{
	  /*printf("Now it will be turned into a picture!\n");*/
	  count++;

	  make_pic_from_region(current, n++, begbit, endbit);
	}
      else
	count++;

      current=current->nregion;
    }

  strcpy(reg_fname, "REGIONS_file");
  if((fp=fopen(reg_fname,"w"))==NULL)
    printf("Couldn't open %s for writing.\n", reg_fname);
  else
    {
      fprintf(fp,"%s",command_string);
      fclose(fp);
    }
}

void make_pic_from_region(struct regionstruct *preg, int n, int begbit, int endbit)
{
  int i,j,k,l,width,height,minx,maxx,miny,maxy;
  struct pixelstruct *pix;

  width = max_width(preg,&minx,&maxx);   /* this is only importantif regions are not */
  height = max_height(preg,&miny,&maxy); /*                                  squares */
  /* printf("minx:%d, maxx:%d, miny:%d, maxy:%d (width:%d height:%d)\n", minx, maxx, miny, maxy, width, height); */

  generate_picture(preg,n,width,minx,maxx,height,miny,maxy, begbit, endbit);
}

void generate_picture(struct regionstruct *preg, int n, int w, int minx, int maxx, int h, int miny, int maxy, int begbit, int endbit)
{
  int Kleinbild[w][h],i,j,k,length, colour1, p1, go_on=1;
  int Realbild[w][h];
  struct pixelstruct *pix;
  char fname[STANDARD_STRINGLENGTH];
  char rname[STANDARD_STRINGLENGTH];
  char png_name[STANDARD_STRINGLENGTH];
  char transform_cmd[STANDARD_STRINGLENGTH];
  FILE *fp;

  for(i=0;i<w;i++)
    for(j=0;j<h;j++)
      Kleinbild[i][j] = 254;

  pix=preg->ipixel;

  while(go_on==1)
    {
      colour1 = pix->c1;
      p1 = 0;
      for(k=begbit;k<=endbit;k++)
         p1 += (int)(floor(colour1/pow(2,k))-2*floor(colour1/pow(2,k+1)));

      Kleinbild[pix->px-minx][pix->py-miny] = p1*255/(endbit-begbit+1);
      Realbild[pix->px-minx][pix->py-miny] = pix->c3;

      if(pix->exit_pixel==NULL)
	go_on=0;
      else
	pix = pix->exit_pixel;
    }

  strcpy(directory, REGIONNAMES);
  strcpy(pic_directory, PICTUREDIR);  
                             /* all the little 'real' box-size pictures */
                                        /* are stored here in pgm format */
  strcpy(png_directory, PNGDIR);  /* ... and here in png-format */

  length=length_region(preg);
  sprintf(fname,"%s/qreg-%s%d.%d-%d-%d-%d-%d.pgm",directory,file_stem,n,preg->ipixel->px,preg->ipixel->py,w,h,length);   /* Name of File for Kleinbild[][]  -  means 'little picture'*/
   sprintf(rname,"%s/reg-%s%d.%d-%d-%d-%d-%d.pgm",pic_directory,file_stem,n,preg->ipixel->px,preg->ipixel->py,w,h,length);   /* Name of File for Realbild[][]   -  means 'real picture '*/
   sprintf(png_name,"%s/reg-%s%d.%d-%d-%d-%d-%d.png",png_directory,file_stem,n,preg->ipixel->px,preg->ipixel->py,w,h,length);   /* Name of File for Realbild[][]   -  means 'real picture '*/

  strcat(command_string,fname);
  strcat(command_string,"\n");

  if((fp=fopen(fname,"w"))==NULL)
    {
      printf("Couldn't open %s for writing !\n",fname);
    }
  else
    {
      /*            fprintf(fp,"P2\n# SEG produced region : Length %d\n%d %d\n255\n", length,w,h);*/

      fprintf(fp,"P2\n%d %d\n255\n", w,h);

      for(j=0;j<h;j++)
	{
	  for(i=0;i<w;i++)
	    fprintf(fp,"%d  ", Kleinbild[i][j]);  
	  fprintf(fp,"\n");
	}

      fclose(fp);
    }

  if((fp=fopen(rname,"w"))==NULL)
    {
      printf("Couldn't open %s for writing !\n",fname);
    }
  else
    {
      /*            fprintf(fp,"P2\n# SEG produced region : Length %d\n%d %d\n255\n", length,w,h);*/

      fprintf(fp,"P2\n%d %d\n255\n", w,h);
      
      for(j=0;j<h;j++)
	{
	  for(i=0;i<w;i++)
	    fprintf(fp,"%d  ", Realbild[i][j]);  
	  fprintf(fp,"\n");
	}

      fclose(fp);
    }

  sprintf(transform_cmd, "convert %s %s",  rname, png_name);
  system(transform_cmd);
}

length_region(struct regionstruct *preg)
{
  int i;
  struct pixelstruct *pix;

  i=1;

  pix=preg->ipixel;
  
  while(pix->exit_pixel!=NULL)
    {
      pix=pix->exit_pixel;
      i++;
    }

  /* printf("Length = %d\n", i); */

  return i;
}

max_width(struct regionstruct *preg, int *minx, int *maxx)
{
  int min=10000,max=0, go_on=1;
  struct pixelstruct *pix;

  pix = preg->ipixel;

  while(go_on==1)
    {
      if(pix->px>max)
	max=pix->px;
      if(pix->px<min)
	min=pix->px;

      if(pix->exit_pixel==NULL)
	go_on=0;
      else
	pix=pix->exit_pixel;
    }

  *minx = min;
  *maxx = max;
  /*  printf("Max_x: %d, Min_x:%d\n", max, min); */

  return (max - min + 1);
}

max_height(struct regionstruct *preg, int *miny, int *maxy)
{
  int min=10000,max=0, go_on=1;
  struct pixelstruct *pix;

  pix = preg->ipixel;

  while(go_on==1)
    {
      if(pix->py>max)
	max=pix->py;	  
      if(pix->py<min)
	min=pix->py;

      if(pix->exit_pixel==NULL)
	go_on=0;
      else
	pix=pix->exit_pixel;
    }

  *miny = min;
  *maxy = max;
  /*  printf("Max_y: %d, Min_y:%d\n", max, min);  */

  return (max - min + 1);
}

static int pcompare(struct regionstruct *r1, struct regionstruct *r2)
{
  if (r1->ipixel->c2 > r2->ipixel->c2)
    return (1);
  if (r1->ipixel->c2 < r2->ipixel->c2)
    return (-1);
  return (0);
}

void join(struct regionstruct *region1, struct regionstruct *region2)
{
  struct regionstruct *creg;
  struct pixelstruct *cur;
  
  cur = region1->ipixel;
  while(cur->c1 != -1)
    cur = cur->npixel;
  cur->npixel = region2->ipixel;

  cur = region2->ipixel;
  while(cur->c1!=-1)
    {
      cur->fpixel=region1->ipixel;
      cur->parent=region1->ipixel->parent;
      cur=cur->npixel;
    }

  if((region2!=firstregion)&&(region2!=lastregion))
    {
      region2->formregion->nregion=region2->nregion;
      region2->nregion->formregion=region2->formregion;
    }
  else if(region2==lastregion)
    {
      region2->formregion->nregion=NULL;      
      lastregion=region2->formregion;
    }
  else
    {
      firstregion = region2->nregion;
      firstregion->formregion=NULL;
    }

  free(region2);

  RN = RN-1;
}
