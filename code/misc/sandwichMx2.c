#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*declare variables*/
  const mwSize *dims;
  mwSize *dimsout;
  mwIndex indx;
  int i, j, k, numdims;
  int numelin, numrows;
  mxClassID classid;
  double *input1r, *input1i, *input2r, *input2i, *output1r, *output1i;
  double a[100][2],  b[2][2],  c[100][2];
  double ai[100][2], bi[2][2], ci[100][2];
  
  /*figure out the classid*/
  classid = mxGetClassID(prhs[0]);
     
  /*check inputs*/
  if (nrhs!=2)
    mexErrMsgTxt("Wrong number of input arguments");
  
  /*associate inputs*/
  input1r = mxGetData(prhs[0]);
  input1i = mxGetImagData(prhs[0]);
  
  input2r = mxGetData(prhs[1]);
  input2i = mxGetImagData(prhs[1]);
  
  /*figure out dimension info and number of elements*/
  dims    = mxGetDimensions(prhs[0]);
  numdims = mxGetNumberOfDimensions(prhs[0]);
  numelin = mxGetNumberOfElements(prhs[0]);
  numrows = dims[0];
  
  /* create the vector which contains the dimensionality of the output */
  dimsout = mxMalloc(numdims * sizeof(mwSize));
  for (i=0; i<numdims; i++) {
      dimsout[i] = dims[i];
  }
  dimsout[1] = numrows;
 
  /*associate output*/
  if (input1i == NULL && input2i == NULL)
  {
    plhs[0]  = mxCreateNumericArray(numdims, dimsout, classid, mxREAL);
    output1r = mxGetData(plhs[0]);
  }
  else
  {
    plhs[0]  = mxCreateNumericArray(numdims, dimsout, classid, mxCOMPLEX);
    output1r = mxGetData(plhs[0]);
    output1i = mxGetImagData(plhs[0]);
  }
  
  /* do the computation*/
  if (input1i == NULL && input2i == NULL)
  {  
    for (i=0; i<numelin/(numrows*2); i++)
    {
      
      b[0][0] = input2r[i*4  ];
      b[1][0] = input2r[i*4+1];
      b[0][1] = input2r[i*4+2];
      b[1][1] = input2r[i*4+3];
      
      for (j=0; j<numrows; j++)
      {
        a[j][0] = input1r[i*(numrows*2)+j        ];
        a[j][1] = input1r[i*(numrows*2)+j+numrows];
        
        c[j][0] = a[j][0]*b[0][0]+a[j][1]*b[1][0];
        c[j][1] = a[j][0]*b[0][1]+a[j][1]*b[1][1];
      }
      
      for (j=0; j<numrows; j++)
      {
        for (k=j; k<numrows; k++)
        {
          output1r[i*(numrows*numrows)+j*numrows+k] = c[k][0]*a[j][0]+c[k][1]*a[j][1];
          output1r[i*(numrows*numrows)+k*numrows+j] = output1r[i*(numrows*numrows)+j*numrows+k];
        }
      }    
    }
    return;
  }
  else if (input1i == NULL)
  {
    for (i=0; i<numelin/(numrows*2); i++)
    {
      b[0][0] = input2r[i*4  ];
      b[1][0] = input2r[i*4+1];
      b[0][1] = input2r[i*4+2];
      b[1][1] = input2r[i*4+3];
      
      bi[0][0] = input2i[i*4  ];
      bi[1][0] = input2i[i*4+1];
      bi[0][1] = input2i[i*4+2];
      bi[1][1] = input2i[i*4+3];
      
      for (j=0; j<numrows; j++)
      {
        a[j][0] = input1r[i*(numrows*2)+j        ];
        a[j][1] = input1r[i*(numrows*2)+j+numrows];
        
        c[j][0] = a[j][0]*b[0][0]+a[j][1]*b[1][0];
        c[j][1] = a[j][0]*b[0][1]+a[j][1]*b[1][1];
        
        ci[j][0] = a[j][0]*bi[0][0]+a[j][1]*bi[1][0];
        ci[j][1] = a[j][0]*bi[0][1]+a[j][1]*bi[1][1];
      }
      
      for (j=0; j<numrows; j++)
      {
        for (k=j; k<numrows; k++)
        {
          output1r[i*(numrows*numrows)+j*numrows+k] = c[k][0]*a[j][0]+c[k][1]*a[j][1];
          output1r[i*(numrows*numrows)+k*numrows+j] = output1r[i*(numrows*numrows)+j*numrows+k];
        
          output1i[i*(numrows*numrows)+j*numrows+k] = ci[k][0]*a[j][0]+ci[k][1]*a[j][1];
          output1i[i*(numrows*numrows)+k*numrows+j] = -output1i[i*(numrows*numrows)+j*numrows+k];
        }
      }   
    }
    return;
  }
    
  else
  {  
    
    for (i=0; i<numelin/(numrows*2); i++)
    {
      b[0][0] = input2r[i*4  ];
      b[1][0] = input2r[i*4+1];
      b[0][1] = input2r[i*4+2];
      b[1][1] = input2r[i*4+3];
      
      bi[0][0] = input2i[i*4  ];
      bi[1][0] = input2i[i*4+1];
      bi[0][1] = input2i[i*4+2];
      bi[1][1] = input2i[i*4+3];
      
      for (j=0; j<numrows; j++)
      {
        a[j][0] = input1r[i*(numrows*2)+j        ];
        a[j][1] = input1r[i*(numrows*2)+j+numrows];
        
        ai[j][0] = input1i[i*(numrows*2)+j        ];
        ai[j][1] = input1i[i*(numrows*2)+j+numrows];
        
        c[j][0] = a[j][0]*b[0][0]+a[j][1]*b[1][0]-ai[j][0]*bi[0][0]-ai[j][1]*bi[1][0];
        c[j][1] = a[j][0]*b[0][1]+a[j][1]*b[1][1]-ai[j][0]*bi[0][1]-ai[j][1]*bi[1][1];
        
        ci[j][0] = a[j][0]*bi[0][0]+a[j][1]*bi[1][0]+ai[j][0]*b[0][0]+ai[j][1]*b[1][0];
        ci[j][1] = a[j][0]*bi[0][1]+a[j][1]*bi[1][1]+ai[j][0]*b[0][1]+ai[j][1]*b[1][1];
      }
      
      for (j=0; j<numrows; j++)
      {
        for (k=j; k<numrows; k++)
        {
          output1r[i*(numrows*numrows)+j*numrows+k] = c[k][0]*a[j][0]+c[k][1]*a[j][1]+ci[k][0]*ai[j][0]+ci[k][1]*ai[j][1];
          output1r[i*(numrows*numrows)+k*numrows+j] = output1r[i*(numrows*numrows)+j*numrows+k];
        
          output1i[i*(numrows*numrows)+j*numrows+k] = ci[k][0]*a[j][0]+ci[k][1]*a[j][1]-c[k][0]*ai[j][0]-c[k][1]*ai[j][1];
          output1i[i*(numrows*numrows)+k*numrows+j] = -output1i[i*(numrows*numrows)+j*numrows+k];
        }
      }   
    }
    return;
    
  }
}

/*in the following c and b are swapped with respect to the above
%a b     e f'   a' c'
%c d     f h    b' d'
%
%a*e+b*f  a*f'+b*h  a' c'
%c*e+d*f  c*f'+d*h  b' d'
%
%a*e*a'+b*f*a'+a*f'*b'+b*h*b' a*e*c'+b*f*c'+a*f'*d'+b*h*d'
%c*e*a'+d*f*a'+c*f'*b'+d*h*b' c*e*c'+d*f*c'+c*f'*d'+d*h*d'
%
%e*abs(a)^2 + f*(a'*b) + f'*(a*b') + h*abs(b)^2    e*a*c'    + f*b*c'   + f'*a*d'   + h*b*d'
%e*a'*c    + f*a'*d   + f'*b'*c   + h*b'*d       e*abs(c)^2 + f*(c'*d) + f'*(c*d') + h*abs(d)^2*/

