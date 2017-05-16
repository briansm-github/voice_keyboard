// Voice Keyboard. By Brian Smith.
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <fcntl.h>
#include <sys/soundcard.h>
#include <alsa/asoundlib.h>
#include <termios.h>


#define MAXCB 10000 // maximum stored frames (words times average word size)

#define WINDOW 280 // analysis window (35ms)

#define STEP 80 // for 100fps at 8K (10ms)
#define ORDER 10 // LPC order
#define FFT 32
#define FF2 17
int bin[FF2]={0,16,8,24,4,20,11,19,2,18,10,21,5,25,9,17,1};

#define EMPH 0.9375 // pre-emphasis level
#define TOL 1.5 // dB above ambient to trigger start/stop
#define MINSIZE 15  // minimum acceptible sample size in frames (x10ms)
#define MAXSIZE 70 // maximum acceptible sample size in frames (x10ms)

// globals...
float realtwiddle[FF2],imtwiddle[FF2];
float cb[MAXCB][FF2];
char cbp[MAXCB];
int start[MAXCB],size[MAXCB],cbsize=0;
int parsed=0,trained=0;

float word[500][FF2],safeword[500][FF2];
int wsize=0,safewsize=0;

float lpc[ORDER+1];

//----------------------------------------------------------------------
// read keyboard, returns -1 if no key, otherwise key code.
int getkey() {
    int character;
    struct termios orig_term_attr;
    struct termios new_term_attr;

    /* set the terminal to raw mode */
    tcgetattr(fileno(stdin), &orig_term_attr);
    memcpy(&new_term_attr, &orig_term_attr, sizeof(struct termios));
    new_term_attr.c_lflag &= ~(ECHO|ICANON);
    new_term_attr.c_cc[VTIME] = 0;
    new_term_attr.c_cc[VMIN] = 0;
    tcsetattr(fileno(stdin), TCSANOW, &new_term_attr);

    /* read a character from the stdin stream without blocking */
    /*   returns EOF (-1) if no character is available */
    character = fgetc(stdin);

    /* restore the original terminal attributes */
    tcsetattr(fileno(stdin), TCSANOW, &orig_term_attr);

    return character;
}

//----------------------------------------------------------------------
void setup_twiddles(float *realtwiddle, float *imtwiddle, int N)
{
    int n;

    for(n=1; n<N/2; n++)
    {
       realtwiddle[n] = cos(2.0*M_PI*n/N);
       imtwiddle[n] = -sin(2.0*M_PI*n/N);
    }
}

//-----------------------------------------------------------------------
void simple_fft(float *real, float *im, float *realtwiddle, float *imtwiddle, int N)
{
    unsigned int even, odd, span, log=0, rootindex;    // indexes
    float temp;
    int i;

    for(span=N>>1; span; span>>=1, log++)   
    {
        for(odd=span; odd<N; odd++)         // iterate over the dual nodes
        {
	    
            odd |= span;                    // iterate over odd blocks only
            even = odd ^ span;              // even part of the dual node pair
                      
            temp = real[even] + real[odd];       
            real[odd] = real[even] - real[odd];
            real[even] = temp;
           
            temp = im[even] + im[odd];           
            im[odd] = im[even] - im[odd];
            im[even] = temp;
           
            rootindex = (even<<log) & (N-1); // find root of unity index
            if(rootindex)                    // skip rootindex[0] (has an identity)
            {
                temp=realtwiddle[rootindex]*real[odd]-imtwiddle[rootindex]*im[odd];
                im[odd]=realtwiddle[rootindex]*im[odd]+imtwiddle[rootindex]*real[odd];
                real[odd] = temp;
            }
   
        } // end of loop over n
     } // end of loop over FFT stages
} //end of function


//----------------------------------------------------------------------------

void autocorrelate(
  float Sn[],	/* frame of Nsam windowed speech samples */
  float Rn[],	/* array of P+1 autocorrelation coefficients */
  int Nsam,	/* number of windowed samples to use */
  int order	/* order of LPC analysis */
)
{
  int i,j;	/* loop variables */
  float x;

  for(j=0; j<order+1; j++) {
    x = 0.0;
    for(i=0; i<Nsam-j; i++) x+=Sn[i]*Sn[i+j];
    Rn[j]=x;
  }
}

//----------------------------------------------------------------------------
// Levinson-Durbin function based on Speex 1.0.5 lpc.c....
void wld(
    float       * lpc, /*      [0...p-1] LPC coefficients      */
    const float * ac,  /*  in: [0...p] autocorrelation values  */
    int p
    )
{
   int i, j;  float r, error = ac[0];

   for (i = 0; i < p; i++) {

      /* Sum up this iteration's reflection coefficient.
       */
      r = -ac[i + 1];
      for (j = 0; j < i; j++) r -= lpc[j] * ac[i - j];
      r /= error;

      /*  Update LPC coefficients and total error.
       */
      lpc[i] = r;
      for (j = 0; j < i/2; j++) {
         float tmp  = lpc[j];
         lpc[j]     += r * lpc[i-1-j];
         lpc[i-1-j] += r * tmp;
      }
      if (i % 2) lpc[j] += lpc[j] * r;

      error *= 1.0 - r * r;
   }
}

//---------------------------------------------------------------
// add a word into the codebook.
int train(int p)
{
  FILE *fp;
  int i,n,k;
  
  bcopy(word,&cb[start[cbsize]][0],wsize*FF2*sizeof(int));
  size[cbsize]=wsize; start[cbsize+1]=start[cbsize]+wsize;
  if (p==10) p='R'; // return key
  if (p==127) p='B'; // backspace key
  if (p==32) p='S'; // space key
  cbp[cbsize]=p;  
  fprintf(stderr,"\n*** TRAINED %i or '%c' size=%i\n",cbsize,p,wsize);
  cbsize++;
  
  // write updated codebook back to disk
  fp=fopen("cb.txt","w");
  for (i=0; i<cbsize; i++) {
    fprintf(fp,"%c\n",cbp[i]);
    fprintf(fp,"%i\n",size[i]);
    for (n=0; n<size[i]; n++) {
      for (k=0; k<FF2; k++) fprintf(fp,"%i,",(int)cb[start[i]+n][k]);
      fprintf(fp,"\n");
    }
  } 
  fclose(fp); 
  trained++;
}

//---------------------------------------------------------------
// search codebook, see what word[] most closely matches...
int closest()
{
   FILE *fpout;
   float dist,bestdist,d;
   int best=0,c,i,n,minsize;
   
   bestdist=9999999999.0;
   for (c=0; c<cbsize; c++) {
     dist=0;
     minsize=wsize; if (size[c]<minsize) minsize=size[c];
     for (i=0; i<minsize; i++) {
       for (n=0; n<FF2; n++)
        {d=word[i][n]-cb[start[c]+i][n]; dist+=d*d;}
     }
     dist/=minsize; // average out the distance over the framecount
     if (dist<bestdist && abs(wsize-size[c])<10) {bestdist=dist; best=c;}
   }
   fprintf(stderr,"BEST=%i	or %c, distance=%.f, wsize=%i cbsize=%i\n",
           best, cbp[best],bestdist,wsize,size[best]);
   // store away safe copy of this known-good word...
   safewsize=wsize; bcopy(word,safeword,wsize*FF2*sizeof(int));
   parsed++;
   // only write output if last output was already dealt with fully...
   if (!(fpout=fopen("/dev/shm/stt","r"))) {
     fpout=fopen("/dev/shm/stt","w");
     fprintf(fpout,"%c",cbp[best]); fclose(fpout);
   } else fclose(fpout);
   return(cbp[best]);
}

//---------------------------------------------------------------

int main()
{ 

  FILE *fp;
  snd_pcm_t *handle;
  snd_pcm_hw_params_t *hw_params;
  int rate=8000;
  float f[WINDOW],hann[WINDOW],w[WINDOW],s0,s1=0,tot;
  float ac[ORDER+1],lcp[ORDER+1];
  short buf[1000],sample,s[WINDOW];
  int i,j,n,k,b,c,best,key;
  
  double e;
  double min;
  int silcount=0;

  int sound=0; // boolean start/stop
  float f2[FFT];
  float real[FFT],imag[FFT];

  float amp[WINDOW],pha[WINDOW];
  int frame=0;
  float bestdist,dist; 
   
  setup_twiddles(realtwiddle,imtwiddle,FFT);
  
  // Set up Hann window...
  for(i=0; i<WINDOW; i++) hann[i]=0.5-0.5*cos(2.0*M_PI*i/(WINDOW-1));
    
  cbsize=0; start[0]=0; size[0]=0; cbp[0]='?';
  if (fp=fopen("cb.txt","r")) {
    while(!feof(fp)) {
      fscanf(fp,"%c\n",&cbp[cbsize]);
      fscanf(fp,"%i\n",&size[cbsize]);
      for (i=start[cbsize]; i<start[cbsize]+size[cbsize]; i++) {
        for (n=0; n<FF2; n++) fscanf(fp,"%f,",&cb[i][n]); fscanf(fp,"\n");
      }
      start[cbsize+1]=start[cbsize]+size[cbsize];  cbsize++; 
    }
    fclose(fp); printf("read codebook, size=%i\n",cbsize);
  }
 
  //---------------------------------

  fp=fopen("/tmp/output.raw","w"); // a dump of recorded audio if required

  // set up microphone via ALSA....
  snd_pcm_open(&handle, "default", SND_PCM_STREAM_CAPTURE, 0);	
  snd_pcm_hw_params_malloc(&hw_params);			 
  snd_pcm_hw_params_any(handle, hw_params);
  snd_pcm_hw_params_set_access(handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED);	
  snd_pcm_hw_params_set_format(handle, hw_params, SND_PCM_FORMAT_S16_LE);	
  snd_pcm_hw_params_set_rate_near(handle, hw_params, &rate, 0);
  snd_pcm_hw_params_set_channels(handle, hw_params, 2);
  snd_pcm_hw_params(handle, hw_params);
  snd_pcm_hw_params_free(hw_params);
  snd_pcm_prepare(handle);

  // allow mic to settle...
  snd_pcm_readi(handle, buf, WINDOW); for (i=0; i<WINDOW; i++) s[i]=buf[i*2];
  s1=s[WINDOW-1];
  // first, fill up the initial window...
  snd_pcm_readi(handle, buf, WINDOW); for (i=0; i<WINDOW; i++) s[i]=buf[i*2];
  for (i=0; i<WINDOW; i++) {
    sample=s[i]; s0=(float)sample; f[i]=s0-s1*EMPH;  s1=s0;
  }
    
  while(1) {

    for (i=0; i<WINDOW-STEP; i++) f[i]=f[i+STEP]; // shift samples down
   
    snd_pcm_readi(handle, buf, STEP); for (i=0; i<STEP; i++) s[i]=buf[i*2];
    // fwrite(s,STEP,2,fp); // write out recorded audio for diagnostics?
   
    for (i=WINDOW-STEP,j=0; i<WINDOW; i++,j++) {
      sample=s[j]; s0=(float)sample;   s0/=32768.0;
      f[i]=s0-s1*EMPH; s1=s0; // pre-emphasis  
    }
    for (i=0; i<WINDOW; i++) w[i]=f[i]; 
    
    // remove any DC level.... (probably only needed for crappy analog mics)
    tot=0; for (i=0; i<WINDOW; i++) tot+=w[i]; tot/=WINDOW;
    for (i=0; i<WINDOW; i++) w[i]-=tot;
    
    for (i=0; i<WINDOW; i++) w[i]*=hann[i]; // apply window to frame
    
    autocorrelate(w,ac,WINDOW,ORDER);
    if (ac[0]!=0.0) {
      e=log(sqrtf(ac[0]));  // energy is just log-RMS
      wld(&lpc[1],ac,ORDER); lpc[0]=1.0;
    }
    if (frame==0) min=e; // initial ambient noise estimate
    //printf("e=%f\n",e);
    
    // get (bounded) log-power 17bin Fourier transform of the LPC10 filter....
    for (i=0; i<FFT; i++) {real[i]=0; imag[i]=0;}
    for (i=0; i<=ORDER; i++) real[i]=lpc[i];
    simple_fft(real,imag,realtwiddle,imtwiddle,FFT);
    for (i=0; i<FF2; i++) {
      b=bin[i];
      f2[i]=1.0/sqrtf(real[b]*real[b]+imag[b]*imag[b]);
      f2[i]=logf(f2[i]);
      if (f2[i]<-1.0) f2[i]=-1.0;
      if (f2[i]>2.0) f2[i]=2.0;
    }
   for (i=1; i<FF2; i++) if (f2[i]>1.5) f2[i]=1.5; // more bounding
    for (i=0; i<FF2; i++) f2[i]*=100.0;

    if (sound==0 && e>min+TOL)  { // start recording...
      sound=1; wsize=0;  
      bcopy(f2,&word[wsize],FF2*4); wsize++;
    }
    else if (sound==1 && e>min+TOL) { // continue reading word...
      bcopy(f2,&word[wsize],FF2*4); wsize++;
    }
    else if (sound==1 && e<min+TOL) silcount++;
    if (sound==1 && silcount==20) { // finished reading word
       //printf("wsize=%i\n",wsize);
       if (wsize>MINSIZE && wsize<MAXSIZE) best=closest();
       else  printf("(rejected %i frames)\n",wsize);
       sound=0;
       min=e;
       silcount=0;
    }

   if (e<min) min=e; // update estimate of ambient noise level
   frame++;  

   key=getkey();
   if (key==27) { // ESCAPE key = exit.
     printf("\nExiting.\n");
     if (parsed!=0) 
       printf("Items parsed=%i, Items trained=%i, Accuracy=%i%\n",
               parsed,trained,(parsed-trained)*100/parsed);
     fclose(fp); exit(0);
   }

   if (key!=-1) {  // restore 'safe' copy of last good word and train it...
      wsize=safewsize; bcopy(safeword,word,wsize*FF2*sizeof(int));
      train(key);
   }
  }
}
