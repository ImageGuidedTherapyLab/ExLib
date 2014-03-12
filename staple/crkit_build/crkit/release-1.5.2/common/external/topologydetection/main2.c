#include "genus0.h"

int main(void)
  {
  int dims[3]={150,100,50};  /* dimensions of test volume */
  int i,totlen;
  unsigned short *input;
  FILE  *fp;  

  genus0parameters g0[1];  /* need an instance of genus0parameters */

  genus0init(g0);  /* initialize the instance, set default parameters */

  /* read in the binary volume, and set g0->dims[0..2] */
  totlen=1;
  for (i=0;i<3;i++) totlen*=(g0->dims[i]=dims[i]); 
  input=(unsigned short *)calloc(totlen,sizeof(unsigned short));  
  fp = fopen ("Torus.bin", "rb" );
  fread(input,sizeof(unsigned short),totlen,fp);
  fclose(fp);
  for (i = 0; i < totlen; i++) {
    input[i] = (input[i] != 0);
  }

  /* set some parameters/options */
  g0->cut_loops=1;
  g0->connectivity=18;
  g0->connected_component=0;
  g0->input=input;
  g0->value=1;
  g0->alt_value=2;
  g0->contour_value=3;
  g0->alt_contour_value=4;
  g0->any_genus=0;
  g0->biggest_component=0;
  g0->pad[0]=g0->pad[1]=g0->pad[2]=2;
  g0->ijk2ras=NULL;
  g0->verbose=1;
  g0->extraijkscale[2]=1;

  /* call the function! */
  if (genus0(g0)) return(1); /* check for error */

  /* write out the results */
  if(g0->verbose) printf("Saving data files...\n");
  fp = fopen ("tris.bin", "wb" );
  fwrite(g0->triangles,sizeof(int),g0->tri_count*3,fp);
  fclose(fp);
  fp = fopen ("verts.bin", "wb" );
  fwrite(g0->vertices,sizeof(float),g0->vert_count*3,fp);
  fclose(fp);
  fp = fopen ("outvol.bin", "wb" );
  fwrite(g0->output,sizeof(unsigned short),totlen,fp);
  fclose(fp);

  genus0destruct(g0);  /* free up what's left in memory  */

  return(0);
  }
