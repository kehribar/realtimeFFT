/**************************************************************************
* Real time FFT + time domain graph with openGL
*
* ihsan Kehribar - March 2013
**************************************************************************/
#include <stdlib.h>
#ifdef __APPLE__
   #include <GLUT/glut.h>
#else
   #include <GL/glut.h>
#endif
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "portaudio.h"
#include "fdacoefs.h" // Coefficients for 44.1kHz sample rate, ~2kHz cutoff low pass
#include "fftw3.h"

#define N 16384 
#define REAL 0
#define IMAG 1
#define F_SAMP 44100
#define SAMPLE_RATE  F_SAMP
#define FRAMES_PER_BUFFER 1024
#define FIRSIZE 101 
#define FILTER_ENABLE 1
#define SUBSAMPLE 1
#define SAMPLECOUNT (FRAMES_PER_BUFFER / SUBSAMPLE)
#define PI 3.14159265358979

/* default window size */
const int HEIGHT = 450;
const int WIDTH = 1000;

PaStreamParameters inputParameters;
PaStream *stream = NULL;
PaError err;
fftw_complex *in, *out, *prein; /* pointers for FFT input, output */
fftw_plan my_plan; /* store the type of FFT we want to perform */
int errCount = 0;
float xcirc[FIRSIZE]; 
uint16_t newest = 0;
float sampleBlock_sub[SAMPLECOUNT];
int i,q,t;
double power;
double amplitude;
double freq_step = (double)F_SAMP / (double)N / SUBSAMPLE;
float *sampleBlock;
int numBytes;    
float hammingWindow[N];
float powArray[N];
char printBuffer[128];

int xCursor1;
int xCursor2;
int yCursor1;
int yCursor2;

float updateFir(float newSample)
{
    uint16_t k;
    float temp1;
    float temp2;
    int32_t x_index;
    float mul_result;
    float accumulator = 0.0f;
    
    /* Circularly increment newest index */
    ++newest;

    /* Fold if neccessary ... */
    if(newest == FIRSIZE)
    {
        newest = 0;
    }

    /* Put new sample in delay line */
    xcirc[newest] = newSample;

    /* Do convolution sum */
    x_index = newest;
    for(k = 0; k < FIRSIZE; k++)
    {
        /* Do the math ... */
        temp1 = B_fir[k];
        temp2 = xcirc[x_index];
        mul_result = temp1 * temp2;
        accumulator += mul_result;  
    
        /* Circularly decrement x_index */
        --x_index;

        /* Fold if neccessary */
        if(x_index == -1)
        {
            x_index = FIRSIZE-1;
        }
    }

    /* Return the filter result ... */
    return accumulator;
}   

void createHamming(float* array,int size)
{
   int i = 0;
   for(i=0;i<size;i++)
   {
      array[i] = 0.54 - (0.46 * cos((2*PI*i)/(size-1)));
   }
}

static void idle_function(void)
{
   err = Pa_ReadStream( stream, sampleBlock, FRAMES_PER_BUFFER );
   
   if(err != paNoError)
      printf("Error count: %d\tProblem: %s\n",++errCount,Pa_GetErrorText(err));

   #if FILTER_ENABLE
      /* apply the low pass filter */
      for(i=0;i<FRAMES_PER_BUFFER;i++)
      {
        sampleBlock[i]= updateFir(sampleBlock[i]);
      }
   #endif

   /* subsampling ... */
   for(i=0;i<SAMPLECOUNT;i++)
   {
     sampleBlock_sub[i] = sampleBlock[SUBSAMPLE*i];
   }

   /* shift the input buffer */
   for(i=0;i<(N-SAMPLECOUNT);i++)
   {
     prein[i][REAL] = prein[i+SAMPLECOUNT][REAL];
     prein[i][IMAG] = 0;
   }

   /* add the new samples */
   for(i=N-SAMPLECOUNT;i<N;i++)
   {
     prein[i][REAL] = sampleBlock_sub[i-N+SAMPLECOUNT];
     prein[i][IMAG] = 0;
     //fprintf(record,"%f\n",sampleBlock_sub[i-N+SAMPLECOUNT]);
   }

   /* apply the window function to the samples */
   for(i=0;i<N;i++)
   {
     in[i][REAL] = prein[i][REAL] * hammingWindow[i];
     in[i][IMAG] = 0;
   }

   /* run the fft */
   fftw_execute(my_plan);

   for(i=0;i<((N/2)+1);i++)
   {
     power = (out[i][0]*out[i][0]) + (out[i][1]*out[i][1]); 
     amplitude = sqrt(power);
     amplitude = amplitude/(N/2);
     amplitude *= 1000; // arbitrary gain
     powArray[i]=amplitude;   
   }
      
   glutPostRedisplay();
}

void glutPrint(float x, float y, char* text, float r, float g, float b, float a)
{ 
    if(!text || !strlen(text)) return; 
    bool blending = false; 
    if(glIsEnabled(GL_BLEND)) blending = true; 
    glEnable(GL_BLEND); 
    glColor4f(r,g,b,a); 
    glRasterPos2f(x,y); 
    while (*text) { 
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *text); 
        text++; 
    } 
    if(!blending) glDisable(GL_BLEND); 
}  

void myMouseFunc( int button, int state, int x, int y ) 
{
   printf("button: %d state: %d x_pos: %d y_pos: %d\n",button,state,x,y);
   if(button == 0)
   {
      xCursor1 = x;
      yCursor1 = y;
   }
   else if(button==2)
   {
      xCursor2 = x;
      yCursor2 = y;
   }
}

void keyb(unsigned char key, int x, int y)
{
   if(key==32)
   {
      xCursor1 = 0;
      xCursor2 = 0;
      yCursor1 = 0;
      yCursor2 = 0;
   }
   else
      printf("%d\n",key);  
}

static void display_function(void)
{

   int x_step = 1;
   int step = (N/2) / WIDTH;
   int gain = 5;
   int maximum = 0;
   int max_index;
   int BASE = 2;

   glClear(GL_COLOR_BUFFER_BIT);
   glMatrixMode(GL_MODELVIEW);       
   glLoadIdentity();       

   sprintf(printBuffer,"Freq step: %.2f Hz",freq_step);
   glutPrint(WIDTH/128, HEIGHT-15, printBuffer, 0.85f, 0.85f, 0.85f, 0.0f);
   sprintf(printBuffer,"F_sampling: %.2f Hz",(float)F_SAMP/SUBSAMPLE);
   glutPrint(WIDTH/128, HEIGHT-30, printBuffer, 0.85f, 0.85f, 0.85f, 0.0f);
   sprintf(printBuffer,"Window witdh: %.2f Sec",(float)N/((float)F_SAMP/(float)SUBSAMPLE));
   glutPrint(WIDTH/128, HEIGHT-45, printBuffer, 0.85f, 0.85f, 0.85f, 0.0f);

   glBegin(GL_LINES);     

      glColor3f ( 0.75, 0.75, 0.75);       
      glVertex2f( 0, BASE);         
      glVertex2f( WIDTH, BASE); 

      glColor3f ( 0.95, 0.85, 0.85);       
      glVertex2f( xCursor1, 0);         
      glVertex2f( xCursor1, HEIGHT/4); 

      glColor3f ( 0.75, 0.85, 0.85);       
      glVertex2f( xCursor2, 0);         
      glVertex2f( xCursor2, HEIGHT/4); 

      for(i=0;i<WIDTH;i+=x_step)
      {
         glColor3f ( 1.0, 0.0, 0.0);       
         glVertex2f( i, BASE + (gain * powArray[i*step]));         
         glVertex2f( i+x_step, BASE + (gain * powArray[(i+x_step)*step])); 

         if(powArray[i*step]>maximum)
         {
            maximum = powArray[i*step];
            max_index = i*step;
         }      
      
         glColor3f ( 0.5, 0.5, 0.5);       
         glVertex2f( i, (HEIGHT/2) + (100 * gain * prein[i*step/2][REAL]));         
         glVertex2f( i+x_step, (HEIGHT/2) + (100 * gain * prein[(i+x_step)*step/2][REAL]));  
      }     

      glColor3f ( 0.85, 0.85, 0.85);       
      glVertex2f( (max_index/step), (gain * powArray[max_index]));         
      glVertex2f( (max_index/step), 20 + (gain * powArray[max_index])); 

   glEnd();

   sprintf(printBuffer,"Amplitude: %f",powArray[max_index]);
   glutPrint(max_index/step, 40 + (gain * powArray[max_index]), printBuffer, 0.05f, 0.05f, 0.05f, 0.0f);
   sprintf(printBuffer,"Peak freq: %.2f Hz",(float)max_index*freq_step);
   glutPrint((max_index/step), 60 + (gain * powArray[max_index]), printBuffer, 0.05f, 0.05f, 0.05f, 0.0f);

   int currsor1_bottom;

   if(xCursor1 != 0)
   {
      sprintf(printBuffer,"Cursor Amplitude: %f",powArray[xCursor1*step]);
      glutPrint(WIDTH/128,HEIGHT-60, printBuffer, 0.95f, 0.85f, 0.85f, 0.0f);
      sprintf(printBuffer,"Cursor freq: %.2f Hz",(float)xCursor1*freq_step*step);
      glutPrint(WIDTH/128,HEIGHT-75, printBuffer, 0.95f, 0.85f, 0.85f, 0.0f);
      currsor1_bottom = 75;   
   }
   else
   {
      sprintf(printBuffer,"Cursor Off");
      glutPrint(WIDTH/128,HEIGHT-60, printBuffer, 0.95f, 0.85f, 0.85f, 0.0f);
      currsor1_bottom = 60;
   }

   if(xCursor2 != 0)
   {
      sprintf(printBuffer,"Cursor Amplitude: %f",powArray[xCursor2*step]);
      glutPrint(WIDTH/128,HEIGHT-(currsor1_bottom+15), printBuffer, 0.75f, 0.85f, 0.85f, 0.0f);
      sprintf(printBuffer,"Cursor freq: %.2f Hz",(float)xCursor2*freq_step*step);
      glutPrint(WIDTH/128,HEIGHT-(currsor1_bottom+30), printBuffer, 0.75f, 0.85f, 0.85f, 0.0f);
   }
   else
   {
      sprintf(printBuffer,"Cursor Off");
      glutPrint(WIDTH/128,HEIGHT-(currsor1_bottom+15), printBuffer, 0.75f, 0.85f, 0.85f, 0.0f);
   }

   glutSwapBuffers();
}

void myinit(void)
{
   glClearColor(0.7, 0.7, 0.7, 0.0); // gray background
   glMatrixMode(GL_PROJECTION);      
   glLoadIdentity();                 
   gluOrtho2D( 0, WIDTH, 0, HEIGHT); // defining the corner points of the window
   glMatrixMode(GL_MODELVIEW);       
}

int main(int argc, char** argv)
{  
   in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N); /* memory allocation */
   prein = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N); /* memory allocation */
   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N); /* pointers in and out are of type fftw_complex*/
   my_plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE); /* create the FFT plan */

   numBytes = FRAMES_PER_BUFFER * sizeof(paFloat32) ;
   sampleBlock = (float *) malloc( numBytes );
   
   err = Pa_Initialize();

   inputParameters.device = Pa_GetDefaultInputDevice();
   inputParameters.channelCount = 1;
   inputParameters.sampleFormat = paFloat32;
   inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency ;
   inputParameters.hostApiSpecificStreamInfo = NULL;

   err = Pa_OpenStream(
           &stream,
           &inputParameters,
           NULL,
           SAMPLE_RATE,
           FRAMES_PER_BUFFER,
           paClipOff, 
           NULL, /* no callback, use blocking API */
           NULL /* no callback, so no callback userData */
   );

   err = Pa_StartStream( stream );

   createHamming(hammingWindow,N);

   for(i=0;i<N;i++)
   {
      prein[i][REAL] = 0;
      prein[i][IMAG] = 0;
   }

   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
   glutInitWindowSize(WIDTH, HEIGHT);
   glutCreateWindow("Real time FFT graph");
   glutIdleFunc(&idle_function);
   glutDisplayFunc(&display_function);
   glutMouseFunc(&myMouseFunc);
   glutKeyboardFunc(keyb);
   glEnable(GL_BLEND);

   myinit();
   glutMainLoop();
   
   return 0;
}

