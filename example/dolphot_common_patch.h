#include <fits.h>
#include <unistd.h>
#include <time.h>

#define DOLPHOT_MAIN
#include "dolphot_defs.h"

// Definitions
typedef char fntype[161];
typedef int int2[2];
typedef struct {
   fntype base;
   dpostype dpos;
   dpostype ref2img;
   apsftype apsf;
   float apdata[3];
   double RAper,RChi,RSky0,RSky1,RSky20,RSky21;
   int RPSF;
} imgdatatype;
typedef struct{
   int flag,type;
   double x,y;
   float s,ss,chi,sh,sky;
   float *is,*iss,*icm,*pa,*pb,*pc;
} stype;

#include "dolphot.h"

// Declarations
#ifdef USEACS
#include "acs/acsphot.h"
#endif
#ifdef USEWFPC2
#include "wfpc2/wfpc2phot.h"
#endif
#ifdef USEWFC3
#include "wfc3/wfc3phot.h"
#endif
#ifdef USEROMAN
#include "roman/romanphot.h"
#endif
#ifdef USENIRCAM
#include "nircam/nircamphot.h"
#endif
#ifdef USENIRISS
#include "niriss/nirissphot.h"
#endif
#ifdef USEMIRI
#include "miri/miriphot.h"
#endif

// Input parameters
imgdatatype imgdata[MAXNIMG];
fntype psfstars="",xytfile="",xytpsf="",FakeStars="",FakeOut="",UsePhot="",DiagPlotType="";
int MaxIT=25,RPSF0=13,RCentroid=1,SecondPass=1,RandomFake=1,FakePad=0,PSFres=1,PSFPhot=1,PSFPhotIt=1,FitSky=1,ApCor=1,Align=1,AlignIter=1,AlignOnly=0,Rotate=1,Force1=0,SkipSky=1,FakeStarPSF=1,SearchMode=1,UseWCS=0,ForceSameMag=0,SubResRef=1,VerboseData=0,NegSky=1;
int XMIN2=0,XMAX2=0,YMIN2=0,YMAX2=0,CUSE=-1,GUSE=-1;
double RAper0=2.5,RChi0=-1.0,SigFind=2.5,SigFinal=3.5,NoiseMult=0.05*0.05,FSat=0.999,PosStep=0.25,dPosMax=3.,RCombine=2.0,SigPSF=10.,psfoff=0.,RSky00=4,RSky10=10,RSky200=-1,RSky210=-1,SkySig=2.25,SigFindMult=0.85,FakeMatch=3.,FakePSF=0.72134752,RBig=0.,AlignTol=0,AlignStep=1;
long MCSEED = 0;

// Derived parameters
int templateexists=0,WARMSTART=0,rMark,rMarkPSF,rpsfMax,rPhotPlusPSFmax;

// Allocated variables
int INTSIZE=4;
int LONGSIZE=4;
int FLOATSIZE=4;
int DOUBLESIZE=8;
int PTRSIZE=4;
stype *stars;
double *iRN,*iAIR;
float *iDMAX,*iDMIN,*iDMIN0,***poff,*sky_val,*sky_unc;
fntype *base;
chiptype *data,*sky,*res;
imtype *skyim,*resim;
double *refmult;
float *refcts,*refmag;
float ***cwt[51][51];
float ***chiwt[51][51];
int2 **circ[51][51];
int *sky_set;

// Parameter parsing routines
void perr(char *str) {
   printf("%s\n",str);
   exit(0);
}

static inline double SQR(double x) {return x*x;}

void initimgdata(void) {
   int i,j,k;

   // default photometry parameters
   Timg=1;
   Nimg=1;
   SubPixel=1;
   FPSF=4;
   EPSF=1;
   MinS=1.0;
   MaxS=9.0;
   MaxE=0.5;
   PSFStep=0.25;
   Zero=25.;
   PSFsol=1;
   WFPC2useCTE=1;
   ACSuseCTE=1;
   WFC3useCTE=1;
   InterpPSFlib=1;
   FlagMask=4;
   CombineChi=0;
   ACSpsfType=0;
   WFC3psfType[0]=WFC3psfType[1]=WFC3psfType[2]=0;
   NIRCAMvega=1;
   NIRISSvega=1;
   MIRIvega=1;
   for (i=0;i<MAXNIMG;i++) {
      imgdata[i].base[0]=0;
      //for (j=0;j<84;j++) imgdata[i].dpos[j]=0.;
      for (j=0;j<40;j++) imgdata[i].dpos[j]=0.;
      imgdata[i].dpos[2]=1.;
      for (j=0;j<40;j++) imgdata[i].ref2img[j]=0.;
      imgdata[i].ref2img[1]=1.;
      imgdata[i].ref2img[12]=1.;
      imgdata[i].ref2img[21]=1.;
      imgdata[i].ref2img[32]=1.;
      for (j=0;j<3;j++) for (k=0;k<6;k++) imgdata[i].apsf[j][k]=0.;
      imgdata[i].apsf[0][0]=imgdata[i].apsf[1][0]=3.;
      imgdata[i].apdata[0]=20.;
      imgdata[i].apdata[1]=30.;
      imgdata[i].apdata[2]=50.;
      if (imgdata[i].apsf[0][0]<=0 || imgdata[i].apsf[1][0]<=0) perr("Initial FWHM must be positive");
      imgdata[i].RAper=-1.0;
      imgdata[i].RChi=-1.0;
      imgdata[i].RSky0=-1.0;
      imgdata[i].RSky1=-1.0;
      imgdata[i].RSky20=-1.0;
      imgdata[i].RSky21=-1.0;
      imgdata[i].RPSF=-1;
   }
   // other globals
   DRIZZLE_BASE=0;
   lastpsftype=-1;
   poffreset=0;
   return;
}

int dolphotparam(char*var,char*val) {
   double x;
   int i,j,k,img=0;
   char *ptr0,*ptr,*ptr2;

   if (!strncasecmp(var,"img",3) && (var[3]=='_' || ((img=strtol(var+3,&ptr0,10))>=0 && ptr0>var+3 && *ptr0=='_'))) {
      if (var[3]=='_') {
	 img=0;
	 i=1;
	 ptr0=var+4;
      }
      else {
	 i=0;
	 ptr0++;
      }
      if (!strcasecmp(ptr0,"file")) {
	 if (img==0 && val[0]) templateexists=1;
	 strcpy(imgdata[img].base,val);
	 if (i) for (i=1;i<MAXNIMG;i++) strcpy(imgdata[i].base,imgdata[0].base);
	 return 1;
      }
      if (!strcasecmp(ptr0,"shift")) {
	 imgdata[img].dpos[0]=strtod(val,&ptr);
	 imgdata[img].dpos[1]=strtod(ptr,&ptr2);
	 if (ptr==ptr2 || *ptr2) perr("img_shift requires two parameters");
	 if (i) for (i=1;i<MAXNIMG;i++) for (j=0;j<2;j++) imgdata[i].dpos[j]=imgdata[0].dpos[j];
	 return 1;
      }
      if (!strcasecmp(ptr0,"xform")) {
	 imgdata[img].dpos[2]=strtod(val,&ptr);
	 imgdata[img].dpos[3]=strtod(ptr,&ptr);
	 imgdata[img].dpos[4]=strtod(ptr,&ptr2);
	 if (ptr==ptr2 || *ptr2) perr("img_xform requires three parameters");
	 if (i) for (i=1;i<MAXNIMG;i++) for (j=2;j<5;j++) imgdata[i].dpos[j]=imgdata[0].dpos[j];
	 return 1;
      }
      if (!strcasecmp(ptr0,"psfa") || !strcasecmp(ptr0,"psfb") || !strcasecmp(ptr0,"psfc")) {
	 if (ptr0[3]=='a' || ptr0[3]=='A') j=0;
	 else if (ptr0[3]=='b' || ptr0[3]=='B') j=1;
	 else j=2;
	 ptr2=val;
	 for (k=0;k<6;k++) {ptr=ptr2; imgdata[img].apsf[j][k]=strtod(ptr,&ptr2);}
	 if (ptr==ptr2 || *ptr2) perr("img_psf requires six parameters");
	 if (i) for (i=1;i<MAXNIMG;i++) for (k=0;k<6;k++) imgdata[i].apsf[j][k]=imgdata[0].apsf[j][k];
	 return 1;
      }
      if (!strcasecmp(ptr0,"ref2img")) {
	 ptr2=val;
	 for (k=0;k<20;k++) {ptr=ptr2; imgdata[img].ref2img[k]=strtod(ptr,&ptr2);}
	 imgdata[img].ref2img[0]=0.0;
	 imgdata[img].ref2img[10]=0.0;
	 if (ptr==ptr2 || *ptr2) perr("ref2img requires 20 parameters");
	 if (i) for (i=1;i<MAXNIMG;i++) for (k=0;k<20;k++) imgdata[i].ref2img[k]=imgdata[0].ref2img[k];
	 return 1;
      }
      if (!strcasecmp(ptr0,"aprad")) {
	 imgdata[img].apdata[0]=strtod(val,&ptr);
	 if (ptr==val || *ptr) perr("img_aprad requires one parameter");
	 if (i) for (i=1;i<MAXNIMG;i++) imgdata[i].apdata[0]=imgdata[0].apdata[0];
	 return 1;
      }
      if (!strcasecmp(ptr0,"apsky")) {
	 imgdata[img].apdata[1]=strtod(val,&ptr);
	 imgdata[img].apdata[2]=strtod(ptr,&ptr2);
	 if (ptr==ptr2 || *ptr2) perr("img_apsky requires two parameters");
	 if (i) for (i=1;i<MAXNIMG;i++) for (j=1;j<3;j++) imgdata[i].apdata[j]=imgdata[0].apdata[j];
	 return 1;
      }
      if (!strcasecmp(ptr0,"RAper")) {
	 if (i) RAper0=strtod(val,&ptr);
	 else imgdata[img].RAper=strtod(val,&ptr);
	 if (ptr==val || *ptr) perr("img_RAper requires one parameter");
	 return 1;
      }
      if (!strcasecmp(ptr0,"RChi")) {
	 if (i) RChi0=strtod(val,&ptr);
	 else imgdata[img].RChi=strtod(val,&ptr);
	 if (ptr==val || *ptr) perr("img_RChi requires one parameter");
	 return 1;
      }
      if (!strcasecmp(ptr0,"RSky")) {
	 if (i) {
	    RSky00=strtod(val,&ptr);
	    RSky10=strtod(ptr,&ptr2);
	 }
	 else {
	    imgdata[img].RSky0=strtod(val,&ptr);
	    imgdata[img].RSky1=strtod(ptr,&ptr2);
	 }
	 if (ptr==ptr2 || *ptr2) perr("img_RSky requires two parameters");
	 return 1;
      }
      if (!strcasecmp(ptr0,"RSky0")) {
	 if (i) RSky00=strtod(val,&ptr);
	 else imgdata[img].RSky0=strtod(val,&ptr);
	 if (ptr==val || *ptr) perr("img_RSky0 requires one parameter");
	 return 1;
      }
      if (!strcasecmp(ptr0,"RSky1")) {
	 if (i) RSky10=strtod(val,&ptr);
	 else imgdata[img].RSky1=strtod(val,&ptr);
	 if (ptr==val || *ptr) perr("img_RSky1 requires one parameter");
	 return 1;
      }
      if (!strcasecmp(ptr0,"RSky2")) {
	 if (i) {
	    RSky200=strtod(val,&ptr);
	    RSky210=strtod(ptr,&ptr2);
	 }
	 else {
	    imgdata[img].RSky20=strtod(val,&ptr);
	    imgdata[img].RSky21=strtod(ptr,&ptr2);
	 }
	 if (ptr==ptr2 || *ptr2) perr("img_RSky2 requires two parameters");
	 return 1;
      }
      if (!strcasecmp(ptr0,"RPSF")) {
	 if (i) RPSF0=strtod(val,&ptr);
	 else imgdata[img].RPSF=strtod(val,&ptr);
	 if (ptr==val || *ptr) perr("img_RPSF requires one parameter");
	 return 1;
      }
   }
   if (!strcasecmp(var,"photsec")) {
      GUSE=strtol(val,&ptr,10);
      if (val==ptr) {
	 GUSE=CUSE=-1;
	 XMIN2=XMAX2=YMIN2=YMAX2=0;
      }
      else {
	 CUSE=strtol(ptr,&ptr,10)-1;
	 XMIN2=strtol(ptr,&ptr,10);
	 YMIN2=strtol(ptr,&ptr,10);
	 XMAX2=strtol(ptr,&ptr,10);
	 YMAX2=strtol(ptr,&ptr,10);
      }
      return 1;
   }
   if (!strcasecmp(var,"xytfile")) {strcpy(xytfile,val); return 1;}
   if (!strcasecmp(var,"xytpsf")) {strcpy(xytpsf,val); return 1;}
   if (!strcasecmp(var,"FPSF")) {
      if (!strcasecmp(val,"gauss")) FPSF=1;
      else if (!strcasecmp(val,"lorentz")) FPSF=2;
      else if (!strcasecmp(val,"lorentz2")) FPSF=3;
      else if (!strcasecmp(val,"g+l")) FPSF=4;
      else perr("FPSF=gauss,lorentz,lorentz2,g+l");
      return 1;
   }
   if (!strcasecmp(var,"psfstars")) {strcpy(psfstars,val); return 1;}
   if (!strcasecmp(var,"FakeStars")) {strcpy(FakeStars,val); return 1;}
   if (!strcasecmp(var,"FakeOut")) {strcpy(FakeOut,val); return 1;}
   if (!strcasecmp(var,"UsePhot")) {strcpy(UsePhot,val); return 1;}
   if (!strcasecmp(var,"DiagPlotType")) {strcpy(DiagPlotType,val); if (DiagPlotType[0] && strcasecmp(DiagPlotType,"gif") && strcasecmp(DiagPlotType,"ps") && strcasecmp(DiagPlotType,"png")) perr("DiagPlotType=PS,GIF,PNG"); return 1;}
   if (!strcasecmp(var,"FakeStar")) {printf("Obsolete parameter: \"FakeStar\"\n"); return 1;}
   i=strtol(val,&ptr,10);
   if (!*ptr) {
      if (!strcasecmp(var,"Nimg")) {
	 Nimg=i;
	 if (Nimg<1 || Nimg>=MAXNIMG) {
	    char str[81];
	    sprintf(str,"1<Nimg<%d",MAXNIMG-1);
	    perr(str);
	 }
	 return 1;
      }
      if (!strcasecmp(var,"PSFPhot")) {PSFPhot=i; if (PSFPhot<0 || PSFPhot>2) perr("PSTPhot=0,1,2"); return 1;}
      if (!strcasecmp(var,"PSFPhotIt")) {PSFPhotIt=i; if (PSFPhotIt<1 || PSFPhotIt>5) perr("PSTPhotIt=1-5"); return 1;}
      if (!strcasecmp(var,"FitSky")) {FitSky=i; if (FitSky<0 || FitSky>4) perr("FitSky=0,1,2,3,4"); return 1;}
      if (!strcasecmp(var,"SkipSky")) {SkipSky=i; if (SkipSky<1) perr("SkipSky>=1"); return 1;}
      if (!strcasecmp(var,"NegSky")) {NegSky=i; if (NegSky<0 || NegSky>1) perr("NegSky=0,1"); return 1;}
      if (!strcasecmp(var,"SecondPass")) {SecondPass=i; if (SecondPass<0) perr("SecondPass>=0"); return 1;}
      if (!strcasecmp(var,"MaxIT")) {MaxIT=i; if (MaxIT<=0) perr("MaxIT>0"); return 1;}
      if (!strcasecmp(var,"ForceSameMag")) {ForceSameMag=i;if (ForceSameMag<0 || ForceSameMag>1) perr("ForceSameMag=0,1"); return 1;}
      if (!strcasecmp(var,"ApCor")) {ApCor=i;if (ApCor<0 || ApCor>1) perr("ApCor=0,1"); return 1;}
      if (!strcasecmp(var,"SearchMode")) {SearchMode=i;if (SearchMode<0 || SearchMode>1) perr("SearchMode=0,1"); return 1;}
      if (!strcasecmp(var,"Force1")) {Force1=i; if (Force1<0 || Force1>1) perr("Force1=0,1"); return 1;}
      if (!strcasecmp(var,"UseWCS")) {UseWCS=i; if (UseWCS<0 || UseWCS>2) perr("UseWCS=0,1,2"); return 1;}
      if (!strcasecmp(var,"Align")) {Align=i; if (Align<0 || Align>4) perr("Align=0-4"); return 1;}
      if (!strcasecmp(var,"AlignIter")) {AlignIter=i; if (AlignIter<1) perr("AlignIter>0"); return 1;}
      if (!strcasecmp(var,"AlignOnly")) {AlignOnly=i; if (AlignOnly<0 || AlignOnly>1) perr("AlignOnly=0,1"); return 1;}
      if (!strcasecmp(var,"Rotate")) {Rotate=i; if (Rotate<0 || Rotate>1) perr("Rotate=0,1"); return 1;}
      if (!strcasecmp(var,"RCentroid")) {RCentroid=i; if (RCentroid<=0) perr("RCentroid>0"); return 1;}
      if (!strcasecmp(var,"RPSF")) {RPSF0=i; if (RPSF0<=0) perr("RPSF>0"); return 1;}
      if (!strcasecmp(var,"SubPixel")) {SubPixel=i; if (SubPixel<=0) perr("SubPixel>0"); return 1;}
      if (!strcasecmp(var,"SubResRef")) {SubResRef=i; if (SubResRef<=0) perr("SubResRef>0"); return 1;}
      if (!strcasecmp(var,"EPSF")) {EPSF=i; if (EPSF<0 || EPSF>1) perr("EPSF=0,1"); return 1;}
      if (!strcasecmp(var,"PSFsol")) {PSFsol=i; if (PSFsol<-1 || PSFsol>2) perr("PSFsol=-1,0,1,2"); return 1;}
      if (!strcasecmp(var,"PSFres")) {PSFres=i; if (PSFres<0 || PSFres>1) perr("PSFres=0,1"); return 1;}
      if (!strcasecmp(var,"FakeStarPSF")) {FakeStarPSF=i; if (FakeStarPSF<0 || FakeStarPSF>1) perr("FakeStarPSF=0,1"); return 1;}
      if (!strcasecmp(var,"RandomSeed")) {MCSEED = i; if (MCSEED<=0) perr("RandomSeed>0"); return 1;}
      if (!strcasecmp(var,"RandomFake")) {RandomFake=i; if (RandomFake<0 || RandomFake>1) perr("RandomFake=0,1"); return 1;}
      if (!strcasecmp(var,"FakePad")) {FakePad=i; return 1;}
      if (!strcasecmp(var,"VerboseData")) {VerboseData=i; if (VerboseData<0 || VerboseData>1) perr("VerboseData=0-1"); return 1;}
      if (!strcasecmp(var,"WFPC2useCTE")) {WFPC2useCTE=i; if (WFPC2useCTE<0 || WFPC2useCTE>1) perr("WFPC2useCTE=0,1"); return 1;}
      if (!strcasecmp(var,"ACSuseCTE")) {ACSuseCTE=i; if (ACSuseCTE<0 || ACSuseCTE>1) perr("ACSuseCTE=0,1"); return 1;}
      if (!strcasecmp(var,"WFC3useCTE")) {WFC3useCTE=i; if (WFC3useCTE<0 || WFC3useCTE>1) perr("WFC3useCTE=0,1"); return 1;}
      if (!strcasecmp(var,"InterpPSFlib")) {InterpPSFlib=i; if (InterpPSFlib<0 || InterpPSFlib>1) perr("InterpPSFlib=0,1"); return 1;}
      if (!strcasecmp(var,"FlagMask")) {FlagMask=i; if (FlagMask<0 || FlagMask>7) perr("FlagMask=0-7"); return 1;}
      if (!strcasecmp(var,"CombineChi")) {CombineChi=i; if (CombineChi<0 || CombineChi>1) perr("CombineChi=0-1"); return 1;}
      if (!strcasecmp(var,"ACSpsfType")) {ACSpsfType=i; if (ACSpsfType<0 || ACSpsfType>1) perr("ACSpsfType=0-1"); return 1;}
      if (!strcasecmp(var,"WFC3IRpsfType")) {WFC3psfType[0]=i; if (WFC3psfType[0]<0 || WFC3psfType[0]>2) perr("WFC3IRpsfType=0-2"); return 1;}
      if (!strcasecmp(var,"WFC3UVISpsfType")) {WFC3psfType[1]=WFC3psfType[2]=i; if (WFC3psfType[1]<0 || WFC3psfType[1]>1) perr("WFC3UVISpsfType=0-1"); return 1;}
      if (!strcasecmp(var,"NIRCAMvega")) {NIRCAMvega=i; if (NIRCAMvega<0 || NIRCAMvega>1) perr("NIRCAMvega=0-1"); return 1;}
      if (!strcasecmp(var,"NIRISSvega")) {NIRISSvega=i; if (NIRISSvega<0 || NIRISSvega>1) perr("NIRISSvega=0-1"); return 1;}
      if (!strcasecmp(var,"MIRIvega")) {MIRIvega=i; if (MIRIvega<0 || MIRIvega>1) perr("MIRIvega=0-1"); return 1;}
   }
   x=strtod(val,&ptr);
   if (!*ptr) {
      if (!strcasecmp(var,"RAper")) {RAper0=x; if (RAper0<=0) perr("RAper>0"); return 1;}
      if (!strcasecmp(var,"RChi")) {RChi0=x; if (RChi0<=0) perr("RChi>0"); return 1;}
      if (!strcasecmp(var,"RSky0")) {RSky00=x; return 1;}
      if (!strcasecmp(var,"RSky1")) {RSky10=x; return 1;}
      if (!strcasecmp(var,"SkySig")) {SkySig=x; if (FitSky==1 && SkySig<1) perr("SkySig>=1"); return 1;}
      if (!strcasecmp(var,"SigFind")) {SigFind=x; if (SigFind<=-999) perr("SigFind>0"); return 1;}
      if (!strcasecmp(var,"SigFindMult")) {SigFindMult=x; if (SigFindMult<=0) perr("SigFindMult>0"); return 1;}
      if (!strcasecmp(var,"SigFinal")) {SigFinal=x; if (SigFinal<=-999) perr("SigFinal>0"); return 1;}
      if (!strcasecmp(var,"NoiseMult")) {NoiseMult=x; if (NoiseMult<0) perr("NoiseMult>=0."); NoiseMult*=NoiseMult; return 1;}
      if (!strcasecmp(var,"FSat")) {FSat=x; if (FSat<=0) perr("FSat>0"); if (FSat>1) FSat=1; return 1;}
      if (!strcasecmp(var,"AlignTol")) {AlignTol=x; if (AlignTol<0) perr("AlignTol>=0"); return 1;}
      if (!strcasecmp(var,"AlignStep")) {AlignStep=x; if (AlignStep<=0) perr("AlignStep>0"); return 1;}
      if (!strcasecmp(var,"Zero")) {Zero=x; return 1;}
      if (!strcasecmp(var,"PosStep")) {PosStep=x; if (PosStep<=0) perr("PosStep>0"); return 1;}
      if (!strcasecmp(var,"dPosMax")) {dPosMax=x; if (dPosMax<=0) perr("dPosMax>0"); return 1;}
      if (!strcasecmp(var,"RCombine")) {RCombine=x; if (RCombine<0) perr("RCombine>=0"); return 1;}
      if (!strcasecmp(var,"SigPSF")) {SigPSF=x; if (SigPSF<0) perr("SigPSF>=0"); return 1;}
      if (!strcasecmp(var,"PSFStep")) {PSFStep=x; if (PSFStep<0) perr("PSFStep>=0"); return 1;}
      if (!strcasecmp(var,"MinS")) {MinS=x; if (MinS<0) perr("MinS>=0"); return 1;}
      if (!strcasecmp(var,"MaxS")) {MaxS=x; if (MaxS<MinS) perr("MaxS>=MinS"); return 1;}
      if (!strcasecmp(var,"MaxE")) {MaxE=x; if (MaxE<0 || MaxE>1) perr("0<=MaxE<=1"); return 1;}
      if (!strcasecmp(var,"psfoff")) {psfoff=x; return 1;}
      if (!strcasecmp(var,"FakeMatch")) {FakeMatch=x; if (FakeMatch<=0.) perr("FakeMatch>0"); return 1;}
      if (!strcasecmp(var,"FakePSF")) {FakePSF=x; if (FakePSF<=0.) perr("FakePSF>0"); FakePSF=SQR(0.4246609*FakePSF); return 1;}
   }
   if (!strcasecmp(var,"RSky")) {
      RSky00=x;
      RSky10=strtod(ptr,&ptr2);
      if (ptr==ptr2 || *ptr2) perr("img_RSky requires two parameters");
      return 1;
   }
   return 0;
}

void alloc_img_data(void) {
   int i,x,y,j,k,img;
   float *ptr;

   if (PSFsol>=0 && PSFStep<=0.) {
      printf("PSF Step size is zero; canceling analytical PSF solution\n");
      PSFsol=-1;
   }
   if (templateexists) {
      Timg=Nimg+1;
      if (*xytfile) {
	 if (*xytpsf) WARMSTART=2;
	 else WARMSTART=1;
      }
   }
   else {
      Timg=Nimg;
      if (*xytfile) WARMSTART=1;
   }
   iGAIN=(double*)calloc(Timg,DOUBLESIZE);
   iEXP=(double*)calloc(Timg,DOUBLESIZE);
   iEXP0=(double*)calloc(Timg,DOUBLESIZE);
   iRN=(double*)calloc(Timg,DOUBLESIZE);
   iDMAX=(float*)calloc(Timg,FLOATSIZE);
   iDMIN=(float*)calloc(Timg,FLOATSIZE);
   iDMIN0=(float*)calloc(Timg,FLOATSIZE);
   iEPOCH=(double*)calloc(Timg,DOUBLESIZE);
   iAIR=(double*)calloc(Timg,DOUBLESIZE);
   apcor=(double*)calloc(Timg,DOUBLESIZE);
   apsize=(double*)calloc(Timg,DOUBLESIZE);
   RAper=(double*)calloc(Timg,DOUBLESIZE);
   RChi=(double*)calloc(Timg,DOUBLESIZE);
   RSky0=(double*)calloc(Timg,DOUBLESIZE);
   RSky1=(double*)calloc(Timg,DOUBLESIZE);
   RSky20=(double*)calloc(Timg,DOUBLESIZE);
   RSky21=(double*)calloc(Timg,DOUBLESIZE);
   RPSF=(int*)calloc(Timg,INTSIZE);
   poff=(float***)calloc(Timg,PTRSIZE);
   rphot=(int*)calloc(Timg,INTSIZE);
   dpos=(dpostype*)calloc(Timg,sizeof(dpostype));
   dpos0=(dpostype*)calloc(Timg,sizeof(dpostype));
   ref2img=(dpostype*)calloc(Timg,sizeof(dpostype));
   hstoffset=(hstoffsettype*)calloc(Timg,sizeof(hstoffsettype));
   wcs=(wcstype*)calloc(Timg,sizeof(wcstype));
   wcsref=(wcsreftype*)calloc(Timg,sizeof(wcsreftype));
   wcsminmax=(wcsreftype*)calloc(Timg,sizeof(wcsreftype));
   wcsminmax_set=0;
   apsky=(apskytype*)calloc(Timg,sizeof(apskytype));
   data=(chiptype*)calloc(Timg,sizeof(chiptype));
   sky=(chiptype*)calloc(Timg,sizeof(chiptype));
   res=(chiptype*)calloc(Timg,sizeof(chiptype));
   dataim=(imtype*)calloc(Timg,sizeof(imtype));
   datahd=(imtype*)calloc(Timg,sizeof(imtype));
   skyim=(imtype*)calloc(Timg,sizeof(imtype));
   resim=(imtype*)calloc(Timg,sizeof(imtype));
   apsf=(apsftype*)calloc(Timg,sizeof(apsftype));
   base=(fntype*)calloc(Timg,sizeof(fntype));
   sky_val=(float*)calloc(Timg,FLOATSIZE);
   sky_unc=(float*)calloc(Timg,FLOATSIZE);
   sky_set=(int*)calloc(Timg,INTSIZE);
   hstmode=(hstmodetype*)calloc(Timg,sizeof(hstmodetype));
   if (!iGAIN || !iEXP || !iEXP0 || !iRN || !iDMAX || !iDMIN || !iDMIN0 || !iEPOCH || !iAIR || !apcor || !apsize || !RAper || !RChi || !RSky0 || !RSky1 || !RSky20 || !RSky21 || !RPSF || !poff || !rphot || !dpos || !dpos0 || !ref2img || !wcs || !wcsref || !wcsminmax || !hstoffset || !apsky || !data || !sky || !res || !dataim || !datahd || !skyim || !resim || !apsf || !base || !sky_set || !sky_val || !sky_unc || !hstmode) merr();
   for (y=0;y<=50;y++) for (x=0;x<=50;x++) {
      circ[y][x]=(int2**)calloc(Timg,PTRSIZE);
      cwt[y][x]=(float***)calloc(Timg,PTRSIZE);
      chiwt[y][x]=(float***)calloc(Timg,PTRSIZE);
      if (!circ[y][x] || !cwt[y][x] || !chiwt[y][x]) merr();
   }
   if (WARMSTART==2) if ((refmult=(double*)calloc(Nimg,DOUBLESIZE))==NULL) merr();
   stars=(stype*)calloc(MAXNSTARS,sizeof(stype));
   ptr=(float*)calloc((size_t)Timg*(size_t)6*(size_t)MAXNSTARS,FLOATSIZE);
   if (!stars || !ptr) merr();
   for (i=0;i<MAXNSTARS;i++) {
      stars[i].is=ptr+(i*6)*Timg;
      stars[i].iss=ptr+(i*6+1)*Timg;
      stars[i].icm=ptr+(i*6+2)*Timg;
      stars[i].pa=ptr+(i*6+3)*Timg;
      stars[i].pb=ptr+(i*6+4)*Timg;
      stars[i].pc=ptr+(i*6+5)*Timg;
   }
   if (WARMSTART==2) {
      refcts=(float*)calloc(MAXNSTARS,FLOATSIZE);
      refmag=(float*)calloc(MAXNSTARS,FLOATSIZE);
      if (!refcts || !refmag) merr();
   }
   for (img=0;img<Timg;img++) {
      hstmode[img].inst=NONE;
      if (img<Nimg) i=img+1;
      else i=0;
      strcpy(base[img],imgdata[i].base);
      for (j=0;j<5;j++) dpos0[img][j]=imgdata[i].dpos[j];
      for (j=0;j<40;j++) ref2img[img][j]=imgdata[i].ref2img[j];
      for (j=0;j<3;j++) for (k=0;k<6;k++) apsf[img][j][k]=imgdata[i].apsf[j][k];
   }
   return;
}

void initimgpars(void) {
   int img,i,j;

   for (img=0;img<Timg;img++) {
      if (img<Nimg) i=img+1;
      else i=0;
      apsize[img]=imgdata[i].apdata[0];
      for (j=0;j<2;j++) apsky[img][j]=imgdata[i].apdata[1+j];
      if (imgdata[i].RAper>0.) RAper[img]=imgdata[i].RAper;
      else RAper[img]=RAper0;
      if (imgdata[i].RChi>0.) RChi[img]=imgdata[i].RChi;
      else if (RChi0>0.) RChi[img]=RChi0;
      else RChi[img]=RAper[img];
      if (imgdata[i].RSky0>0.) RSky0[img]=imgdata[i].RSky0;
      else RSky0[img]=RSky00;
      if (imgdata[i].RSky1>0.) RSky1[img]=imgdata[i].RSky1;
      else RSky1[img]=RSky10;
      if (imgdata[i].RPSF>0) RPSF[img]=imgdata[i].RPSF;
      else RPSF[img]=RPSF0;
      if (imgdata[i].RSky20>0.) RSky20[img]=imgdata[i].RSky20;
      else if (RSky200>0.) RSky20[img]=RSky200;
      else RSky20[img]=RAper[img]+1;
      if (imgdata[i].RSky21>0.) RSky21[img]=imgdata[i].RSky21;
      else if (RSky210>0.) RSky21[img]=RSky210;
      else RSky21[img]=RPSF[img];
      // checks for sky radii
      if (FitSky==1 || FitSky==2) {
	 if (RSky0[img]<RAper[img]+0.5) {
	    RSky0[img]=RAper[img]+0.5;
	    printf("Increasing RSky0 to %f for image %d\n",RSky0[img],i);
	    imgdata[i].RSky0=RSky0[img];
	 }
	 if (RSky1[img]<RSky0[img]+1) {
	    RSky1[img]=RSky0[img]+1;
	    printf("Increasing RSky1 to %f for image %d\n",RSky1[img],i);
	    imgdata[i].RSky1=RSky1[img];
	 }
      }
      if (FitSky==2) {
	 if (RPSF[img]<RAper[img]+2) perr("Error: RPSF needs to be at least RAper+2 if using FitSky=2");
	 if (RSky20[img]<RAper[img]+1) {
	    RSky20[img]=RAper[img]+1;
	    printf("Increasing RSky20 to %f for image %d\n",RSky20[img],i);
	    imgdata[i].RSky20=RSky20[img];
	 }
	 else if (RSky20[img]>RPSF[img]-1) {
	    RSky20[img]=RPSF[img]-1;
	    printf("Decreasing RSky20 to %f for image %d\n",RSky20[img],i);
	    imgdata[i].RSky20=RSky20[img];
	 }
	 if (RSky21[img]>RPSF[img]) {
	    RSky21[img]=RPSF[img];
	    printf("Decreasing RSky21 to %f for image %d\n",RSky21[img],i);
	    imgdata[i].RSky21=RSky21[img];
	 }
	 else if (RSky21[img]<RSky20[img]+1) {
	    RSky21[img]=RSky20[img]+1;
	    printf("Increasing RSky21 to %f for image %d\n",RSky21[img],i);
	    imgdata[i].RSky21=RSky21[img];
	 }
      }
      if (RChi[img]>RAper[img]) {
	 RChi[img]=RAper[img];
	 printf("Reducing RChi to %f for image %d\n",RChi[img],i);
	 imgdata[i].RChi=RChi[img];
      }
      hstoffset[img][0] = atoi(getcardval(dataim+img,"DOL_OFFX",0));
      hstoffset[img][1] = atoi(getcardval(dataim+img,"DOL_OFFY",0));
   }
   return;
}

void initcirc(void) {
   int img,x,y,px,py;
   float fx,fy;

   if (PSFStep<=0.) {
      RBig=RAper0;
      if (RBig<MaxS+1.) RBig=MaxS+1.;
   }
   rpsfMax=0;
   rPhotPlusPSFmax=0;
   for (img=0;img<Timg;img++) {
      rphot[img]=(int)(RAper[img]+0.9999);
      if (RPSF[img]>rpsfMax) rpsfMax=RPSF[img];
      if (rphot[img]+RPSF[img]>rPhotPlusPSFmax) rPhotPlusPSFmax=rphot[img]+RPSF[img];
      poff[img]=(float**)calloc(2*RPSF[img]+1,PTRSIZE);
      if (!poff[img]) merr();
      poff[img]+=RPSF[img];
      for (y=-RPSF[img];y<=RPSF[img];y++) {
	 poff[img][y]=(float*)calloc(2*RPSF[img]+1,FLOATSIZE);
	 if (!poff[img][y]) merr();
	 poff[img][y]+=RPSF[img];
      }
   }
   psf=(float**)calloc(rpsfMax*2+1,PTRSIZE);
   if (!psf) merr();
   psf+=rpsfMax;
   for (y=-rpsfMax;y<=rpsfMax;y++) {
      psf[y]=(float*)calloc(rpsfMax*2+1,FLOATSIZE);
      if (!psf[y]) merr();
      psf[y]+=rpsfMax;
   }
   for (img=0;img<Timg;img++) for (y=0;y<=50;y++) for (x=0;x<=50;x++) {
      fx=(x-25)*0.02;
      fy=(y-25)*0.02;
      circ[y][x][img]=(int2*)calloc(rphot[img]*2+1,2*INTSIZE);
      cwt[y][x][img]=(float**)calloc(rphot[img]*2+1,PTRSIZE);
      chiwt[y][x][img]=(float**)calloc(rphot[img]*2+1,PTRSIZE);
      if (!circ[y][x][img] || !cwt[y][x][img] || !chiwt[y][x][img]) merr();
      circ[y][x][img]+=rphot[img];
      cwt[y][x][img]+=rphot[img];
      chiwt[y][x][img]+=rphot[img];
      for (py=-rphot[img];py<=rphot[img];py++) {
	 cwt[y][x][img][py]=(float*)calloc(rphot[img]*2+1,FLOATSIZE);
	 chiwt[y][x][img][py]=(float*)calloc(rphot[img]*2+1,FLOATSIZE);
	 if (!cwt[y][x][img][py] || !chiwt[y][x][img][py]) merr();
	 cwt[y][x][img][py]+=rphot[img];
	 chiwt[y][x][img][py]+=rphot[img];
	 circ[y][x][img][py][0]=500;
	 circ[y][x][img][py][1]=-500;
	 for (px=-rphot[img];px<=rphot[img];px++) {
	    cwt[y][x][img][py][px]=RAper[img]+.465-sqrt((fx-px)*(fx-px)+(fy-py)*(fy-py));
	    if (cwt[y][x][img][py][px]<0) cwt[y][x][img][py][px]=0;
	    else if (cwt[y][x][img][py][px]>1) cwt[y][x][img][py][px]=1;
	    if (cwt[y][x][img][py][px]>0 && px<circ[y][x][img][py][0]) circ[y][x][img][py][0]=px;
	    if (cwt[y][x][img][py][px]>0 && px>circ[y][x][img][py][1]) circ[y][x][img][py][1]=px;
	    chiwt[y][x][img][py][px]=RChi[img]+.465-sqrt((fx-px)*(fx-px)+(fy-py)*(fy-py));
	    if (chiwt[y][x][img][py][px]<0) chiwt[y][x][img][py][px]=0;
	    else if (chiwt[y][x][img][py][px]>1) chiwt[y][x][img][py][px]=1;
	 }
      }
   }
   return;
}

void freecirc(void) {
   int img,x,y,py;

   for (img=0;img<Timg;img++) {
      for (y=-RPSF[img];y<=RPSF[img];y++) free(poff[img][y]-RPSF[img]);
      free(poff[img]-RPSF[img]);
   }
   for (y=-rpsfMax;y<=rpsfMax;y++) free(psf[y]-rpsfMax);
   free(psf-rpsfMax);
   for (img=0;img<Timg;img++) for (y=0;y<=50;y++) for (x=0;x<=50;x++) {
      for (py=-rphot[img];py<=rphot[img];py++) {
	 free(cwt[y][x][img][py]-rphot[img]);
	 free(chiwt[y][x][img][py]-rphot[img]);
      }
      free(cwt[y][x][img]-rphot[img]);
      free(chiwt[y][x][img]-rphot[img]);
      free(circ[y][x][img]-rphot[img]);
   }
}

// Common utilities
typedef struct {
   long lastoffset;
   FILE *f;
   char fname[801];
   char mode[3];
} reopenableFile;

void fopenagain(reopenableFile*ptr) {
   ptr->f=fopen(ptr->fname,ptr->mode);
   if (ptr->f==NULL) {printf("Error: cannot re-open %s\n",ptr->fname); exit(-1);}
   if (ptr->lastoffset!=0) fseek(ptr->f,ptr->lastoffset,SEEK_SET);
}

void fopenfirst(reopenableFile*ptr,char*fname,char*mode,int open) {
   strcpy(ptr->fname,fname);
   strcpy(ptr->mode,mode);
   ptr->lastoffset=0;
   ptr->f=0;
   if (open==1) fopenagain(ptr);
}

void freclose(reopenableFile*ptr,int check,FILE *fwarn) {
   if (check && feof(ptr->f)) {
      printf("**Error: read past end of %s.  File appears truncated.\n",ptr->fname);
      fflush(stdout);
      if (fwarn) {
	 fprintf(fwarn,"**Error: read past end of %s.  File appears truncated.\n",ptr->fname);
	 fflush(fwarn);
      }
   }
   ptr->lastoffset=ftell(ptr->f);
   fclose(ptr->f);
   ptr->f=0;
}

/*
FILE*fopenagain(int img,char*str,char*type,long off) {
   static FILE *f;
   char fn[161];
   if (img>0) sprintf(fn,"%s%s",base[img-1],str);
   else {printf("Stupid error in fopenagain (%d/%d)\n",img+1,Nimg); exit(-1);}
   f=fopen(fn,type);
   if (f==NULL) {printf("Error: cannot re-open %s\n",fn); exit(-1);}
   if (off!=0) fseek(f,off,SEEK_SET);
   return f;
}

void freclose(FILE*f,long*off) {
   *off=ftell(f);
   fclose(f);
   return;
}
*/

static inline int posOK(int i,int x,int y) {
   return (x>=0 && x<dataim[i].X && y>=0 && y<dataim[i].Y);
}

static inline int dataOK(int i,int x,int y) {
   return (data[i][y][x]>iDMIN[i] && data[i][y][x]<iDMAX[i]);
}

static inline int datafOK(int i,int x,int y) {
   return (data[i][y][x]>iDMIN[i] && data[i][y][x]<iDMAX[i]*FSat);
}

static inline int ppixOK(int i,int x,int y) {
   return (x>=0 && x<dataim[i].X && y>=0 && y<dataim[i].Y && data[i][y][x]>iDMIN[i] && data[i][y][x]<iDMAX[i]);
}

static inline int ppixfOK(int i,int x,int y) {
   return (x>=0 && x<dataim[i].X && y>=0 && y<dataim[i].Y && data[i][y][x]>iDMIN[i] && data[i][y][x]<iDMAX[i]*FSat);
}

void readinfo(FILE *finfo,int ext,int z) {
   char sstr[321],str[321];
   int i,j,k;
   void setwcsminmax(void);

   sprintf(sstr,"EXTENSION %d CHIP %d\n",ext,z+1);
   while (fgets(str,321,finfo) && strcmp(str,sstr));
   if (strcmp(str,sstr)) {
      printf("Cannot locate information from initial photometry .info file\n");
      exit(-1);
   }
   fgets(str,321,finfo);
   if (strcmp(str,"Limits\n")) {printf("Illegal format for .info file; cannot locate limits\n"); exit(-1);}
   fgets(str,321,finfo); // gives photsec data
   while (fgets(str,321,finfo) && *str=='*'); // skip any HST stuff
   if (strcmp(str,"Alignment\n")) {printf("Illegal format for .info file; cannot locate alignment\n"); exit(-1);}
   for (i=0;i<Nimg;i++) {
      if (Align<4) {for (j=0;j<5;j++) fscanf(finfo,"%f",dpos[i]+j);}
      else {for (j=0;j<40;j++) fscanf(finfo,"%f",dpos[i]+j);}
      if (UseWCS==2) {
	 for (j=0;j<4;j++) fscanf(finfo,"%lf",wcsref[i]+j);
	 for (j=0;j<80;j++) fscanf(finfo,"%f",wcs[i]+j);
      }
   }
   fgets(str,321,finfo); // end of last line of alignment data
   fgets(str,321,finfo); // "PSF solution" or "Aperture corrections"
   if (!strcmp(str,"PSF solution\n")) {
      for (i=0;i<Nimg;i++) for (j=0;j<3;j++) for (k=0;k<6;k++) fscanf(finfo,"%lf",apsf[i][j]+k);
      fgets(str,321,finfo); // end of last line of PSF data
      fgets(str,321,finfo); // "Aperture corrections"
   }
   if (strcmp(str,"Aperture corrections\n")) {printf("Illegal format for .info file; cannot locate aperture corrections\n"); exit(-1);}
   for (i=0;i<Nimg;i++) {
      fscanf(finfo,"%lf",apcor+i);
      apcor[i]=pow(10,-0.4*apcor[i]);
   }
   fgets(str,321,finfo); // end of last line of aperture data
   if (UseWCS==2) setwcsminmax();
   return;
}

void read_cardvals(int img) {
   if (!isimage(dataim)) return;
   parsecards(datahd+img,iGAIN+img,iRN+img,iEXP+img,iDMIN+img,iDMAX+img,iEPOCH+img,iAIR+img,iEXP0+img,0,1);
   parsecards(dataim+img,iGAIN+img,iRN+img,iEXP+img,iDMIN+img,iDMAX+img,iEPOCH+img,iAIR+img,iEXP0+img,1,0);
   if (iEXP[img]<=0) {
      fprintf(stderr,"ERROR: image %d has zero or negative exposure time\n",img+1);
      exit(-1);
   }
   iRN[img]*=iRN[img]/iGAIN[img]/iGAIN[img];
   iDMIN0[img]=iDMIN[img];
   if (iDMAX[img]>0. && iDMIN0[img]>-iDMAX[img]) iDMIN[img]=-iDMAX[img];
   return;
}

// PSF routines
int X=0,Y=0;
int POSPSF=1;
double *psfkernel;

void setpsfkernel(void) {
   int i;
   double m;

   psfkernel=(double*)calloc(100001,DOUBLESIZE);
   if (!psfkernel) merr();
   m=1./(SubPixel*SubPixel*M_PI);
   if (FPSF==1) m*=0.693;
   else if (FPSF==3) m*=sqrt(2.)-1.;
   for (i=0;i<=100000;i++) {
      if (FPSF==1) psfkernel[i]=m*exp(-0.01*i*0.693);
      else if (FPSF==2) psfkernel[i]=m*1./(1+0.01*i);
      else if (FPSF==3) psfkernel[i]=m*1./(1+0.01*i*0.41421356)/(1+0.01*i*0.41421356);
      else if (FPSF==4) psfkernel[i]=m*(exp(-0.01*i*0.693)+1./(1+0.01*i));
   }
   return;
}

static inline double evalpsf(double x) {
   int i;
   if (x<0) x=-x;
   x*=100.;
   i=(int)x;
   if (i>=100000) {
      x/=100000.;
      return psfkernel[100000]/x/x;
   }
   return (1+i-x)*psfkernel[i]+(x-i)*psfkernel[i+1];
}

static inline void add_poff(int img,int r) {
   int i,j;
   if (POSPSF) {
      //double a=0,b=0;
      for (j=-r;j<=r;j++) for (i=-r;i<=r;i++) {
	 psf[j][i]+=poff[img][j][i];
	 //a+=psf[j][i];
	 if (psf[j][i]<0.) psf[j][i]=0.;
	 //b+=psf[j][i];
      }
      /*
      if (a!=b) {
	 a/=b;
	 for (j=-r;j<=r;j++) for (i=-r;i<=r;i++) psf[j][i]*=a;
      }
      */
   }
   else {
      for (j=-r;j<=r;j++) for (i=-r;i<=r;i++) psf[j][i]+=poff[img][j][i];
   }
   return;
}

void calc1psf(int img,float x,float y,int r,float sa,float sb,float c,int type,int sub) {
   int i,j,ii,jj,sp1=0;
   static int first=1,first_list=1,lastimg,lastr,lasttype,lastsub,lastHSTmode;
   static float *list,lastx,lasty,lasta,lastb,lastc;
   double a,b,m,dx,dy,sps=1.,dy0;
   int solved=0;

   if (sub && FitSky==2 && r<RPSF[img]) r=RPSF[img];
#ifdef USEWFPC2
   if (hstmode[img].inst==WFPC2 && (Force1 || (sa>MinS && sa<MaxS))) {
      if (!calcwfpc2psf(img,x+hstoffset[img][0],y+hstoffset[img][1],r,(hstmode[img].inst!=lastHSTmode || sub!=lastsub))) return;
      lastsub = sub;
      lastHSTmode = hstmode[img].inst;
      solved=1;
   }
#endif
#ifdef USEACS
   if (hstmode[img].inst==ACS && (Force1 || (sa>MinS && sa<MaxS))) {
      if (!calcacspsf(img,x+hstoffset[img][0],y+hstoffset[img][1],r,(hstmode[img].inst!=lastHSTmode || sub!=lastsub))) return;
      lastsub = sub;
      lastHSTmode = hstmode[img].inst;
      solved=1;
   }
#endif
#ifdef USEWFC3
   if (hstmode[img].inst==WFC3 && (Force1 || (sa>MinS && sa<MaxS))) {
      if (!calcwfc3psf(img,x+hstoffset[img][0],y+hstoffset[img][1],r,(hstmode[img].inst!=lastHSTmode || sub!=lastsub))) return;
      lastsub = sub;
      lastHSTmode = hstmode[img].inst;
      solved=1;
   }
#endif
#ifdef USEROMAN
   if (hstmode[img].inst==ROMAN && (Force1 || (sa>MinS && sa<MaxS))) {
      if (!calcromanpsf(img,x+hstoffset[img][0],y+hstoffset[img][1],r,(hstmode[img].inst!=lastHSTmode || sub!=lastsub))) return;
      lastsub = sub;
      lastHSTmode = hstmode[img].inst;
      solved=1;
   }
#endif
#ifdef USENIRCAM
   if (hstmode[img].inst==NIRCAM && (Force1 || (sa>MinS && sa<MaxS))) {
      if (!calcnircampsf(img,x+hstoffset[img][0],y+hstoffset[img][1],r,(hstmode[img].inst!=lastHSTmode || sub!=lastsub))) return;
      lastsub = sub;
      lastHSTmode = hstmode[img].inst;
      solved=1;
   }
#endif
#ifdef USENIRISS
   if (hstmode[img].inst==NIRISS && (Force1 || (sa>MinS && sa<MaxS))) {
      if (!calcnirisspsf(img,x+hstoffset[img][0],y+hstoffset[img][1],r,(hstmode[img].inst!=lastHSTmode || sub!=lastsub))) return;
      lastsub = sub;
      lastHSTmode = hstmode[img].inst;
      solved=1;
   }
#endif
#ifdef USEMIRI
   if (hstmode[img].inst==MIRI && (Force1 || (sa>MinS && sa<MaxS))) {
      if (!calcmiripsf(img,x+hstoffset[img][0],y+hstoffset[img][1],r,(hstmode[img].inst!=lastHSTmode || sub!=lastsub))) return;
      lastsub = sub;
      lastHSTmode = hstmode[img].inst;
      solved=1;
   }
#endif
   if (!solved) {
      //printf("Computing analaytic PSF (%f %f %f)\n",sa,sb,sb);
      if (SubPixel!=1) {
	 sps=1./SubPixel;
	 sp1=SubPixel-1;
      }
      x-=(int)x+0.5;
      y-=(int)y+0.5;
      if (!first && lastpsftype==0 && img==lastimg && x==lastx && y==lasty && r<=lastr && sa==lasta && sb==lastb && c==lastc && type==lasttype && sub==lastsub && !poffreset && NONE==lastHSTmode) return;
      if (first) {
	 setpsfkernel();
	 first = 0;
      }
      lastpsftype=0;
      lastimg=img;
      lastx=x;
      lasty=y;
      lasta=sa;
      lastb=sb;
      lastc=c;
      lasttype=type;
      lastsub=sub;
      lastr=r;
      lastHSTmode=NONE;
      if (sa<=0 || sb<=0) {
	 printf("Error in PSF generation: negative radius\n");
	 //printf("%f\n",data[-100][-10000][-1000000]);
	 exit(-1);
      }
      if (c>2 || c<-2) {
	 printf("Error in PSF generation: illegal cross term\n");
	 exit(-1);
      }
      c/=0.25*sa*sb;
      a=4./(sa*sa);
      b=4./(sb*sb);
      m=1.;
      m=sqrt(a*b-c*c*0.25);
      if (FPSF==2) m=m/(log(1+RPSF[img]*RPSF[img]*1.5*m));
      else if (FPSF==4) m=m/(1./0.693+log(1+RPSF[img]*RPSF[img]*1.5*m));
      if (SubPixel==1) for (j=-r;j<=r;j++) for (i=-r;i<=r;i++) {
	 dx=x-i;
	 dy=y-j;
	 psf[j][i]=m*evalpsf(a*dx*dx+b*dy*dy+c*dx*dy);
      }
      else if (SubPixel==2) for (j=-r;j<=r;j++) for (i=-r;i<=r;i++) {
         double dxp,dxm,dyp,dym;
	 dxm=x-i-0.5*(1.-sps); dxp=dxm+sps;
	 dym=y-j-0.5*(1.-sps); dyp=dym+sps;
	 psf[j][i]=m*(
	    evalpsf(a*dxm*dxm+b*dym*dym+c*dxm*dym)+evalpsf(a*dxp*dxp+b*dym*dym+c*dxp*dym)
	    +evalpsf(a*dxm*dxm+b*dyp*dyp+c*dxm*dyp)+evalpsf(a*dxp*dxp+b*dyp*dyp+c*dxp*dyp)
	    );
      }
      else if (SubPixel==3) for (j=-r;j<=r;j++) for (i=-r;i<=r;i++) {
         double dxp,dxm,dyp,dym;
	 dx=x-i; dxp=dx+sps; dxm=dx-sps;
	 dy=y-j; dyp=dy+sps; dym=dy-sps;
	 psf[j][i]=m*(
	    evalpsf(a*dxm*dxm+b*dym*dym+c*dxm*dym)+evalpsf(a*dx*dx+b*dym*dym+c*dx*dym)+evalpsf(a*dxp*dxp+b*dym*dym+c*dxp*dym)
	    +evalpsf(a*dxm*dxm+b*dy*dy+c*dxm*dy)+evalpsf(a*dx*dx+b*dy*dy+c*dx*dy)+evalpsf(a*dxp*dxp+b*dy*dy+c*dxp*dy)
	    +evalpsf(a*dxm*dxm+b*dyp*dyp+c*dxm*dyp)+evalpsf(a*dx*dx+b*dyp*dyp+c*dx*dyp)+evalpsf(a*dxp*dxp+b*dyp*dyp+c*dxp*dyp)
	    );
      }
      else for (j=-r;j<=r;j++) for (i=-r;i<=r;i++) {
         psf[j][i]=0;
	 dx=x-i-0.5*(1.-sps);
	 dy=dy0=y-j-0.5*(1.-sps);
	 psf[j][i]+=evalpsf(a*dx*dx+b*dy*dy+c*dx*dy);
	 for (jj=0;jj<sp1;jj++) {
	    dy+=sps;
	    psf[j][i]+=evalpsf(a*dx*dx+b*dy*dy+c*dx*dy);
	 }
	 for (ii=0;ii<sp1;ii++) {
	    dx+=sps;
	    dy=dy0;
	    psf[j][i]+=evalpsf(a*dx*dx+b*dy*dy+c*dx*dy);
	    for (jj=0;jj<sp1;jj++) {
	       dy+=sps;
	       psf[j][i]+=evalpsf(a*dx*dx+b*dy*dy+c*dx*dy);
	    }
	 }
	 psf[j][i]*=m;
      }
   }
   if (type<3) add_poff(img,r);
   if (sub && FitSky==2) {
      int iy,ix,x1,y1,i,N;
      float r2,av;

      if (first_list) {
	 list=(float*)calloc((2*rpsfMax+1)*(2*rpsfMax+1),FLOATSIZE);
	 if (!list) merr();
	 first_list=0;
      }
      N=0;
      iy=(int)y;
      ix=(int)x;
      for (y1=iy-RSky21[img];y1<=iy+RSky21[img];y1++) for (x1=ix-RSky21[img];x1<=ix+RSky21[img];x1++) if (sub==2 || ppixOK(img,x1,y1)) {
	 r2=sqrt((x1-ix)*(x1-ix)+(y1-iy)*(y1-iy));
	 if (r2>=RSky20[img] && r2<=RSky21[img]) list[N++]=psf[y1-iy][x1-ix];
	 //if (img==1) if (r2>=RSky20[img] && r2<=RSky21[img]) printf("%d %d %f %f\n",x1-ix,y1-iy,r2,psf[y1-iy][x1-ix]);
      }
      av=0;
      for (i=0;i<N;i++) av+=list[i];
      if (N>0) av/=N;
      for (x1=-r;x1<=r;x1++) for (y1=-r;y1<=r;y1++) psf[y1][x1]-=av;
      /*
      printf("%f\n",av);
      if (img==1) {
	 for (i=-3;i<=3;i++) {
	    for (j=-3;j<=3;j++) printf("%f ",psf[j][i]);
	    printf("\n");
	 }
	 fflush(stdout);
	 exit(0);
      }
      */
   }
   return;
}

void getpsfpars(int IMG,float x,float y,float*a,float*b,float*c) {
   int img;

   x-=X*0.5;
   y-=Y*0.5;
   for (img=0;img<Timg && (img<Nimg || IMG>=0);img++) if (IMG<0 || img==IMG) {
      a[img]=apsf[img][0][0]+x*apsf[img][0][1]+y*apsf[img][0][2]+x*x*apsf[img][0][3]+y*y*apsf[img][0][4]+x*y*apsf[img][0][5];
      b[img]=apsf[img][1][0]+x*apsf[img][1][1]+y*apsf[img][1][2]+x*x*apsf[img][1][3]+y*y*apsf[img][1][4]+x*y*apsf[img][1][5];
      c[img]=apsf[img][2][0]+x*apsf[img][2][1]+y*apsf[img][2][2]+x*x*apsf[img][2][3]+y*y*apsf[img][2][4]+x*y*apsf[img][2][5];
   }
   return;
}

// used in fake star generation
void shift(int img,double x0,double y0,double *x,double *y,int dir) {
   double r2;
   static double *c,*s,*d2;
   static dpostype*ldpos;
   static int first=1,sdpos,*ifirst;

   if (first) {
      int i;
      sdpos=sizeof(dpostype);
      ldpos=(dpostype*)calloc(Timg,sdpos);
      c=(double*)calloc(Timg,DOUBLESIZE);
      s=(double*)calloc(Timg,DOUBLESIZE);
      d2=(double*)calloc(Timg,DOUBLESIZE);
      ifirst=(int*)calloc(Timg,INTSIZE);
      if (!ldpos || !c || !s || !d2 || !ifirst) merr();
      for (i=0;i<Timg;i++) ifirst[i]=1;
      first=0;
   }
   if (Align<4 && (ifirst[img] || dpos[img][3]!=ldpos[img][3] || dpos[img][4]!=ldpos[img][4])) {
      d2[img]=4.*dpos[img][3]/3./(X*X+Y*Y);
      c[img]=cos(M_PI/180.*dpos[img][4]);
      s[img]=sin(M_PI/180.*dpos[img][4]);
      memcpy(ldpos[img],dpos[img],sdpos);
      ifirst[img]=0;
   }

   // forward transformation (image to reference)
   if (dir<0) {
      if (UseWCS==2) {
	 // move origin to reference pixel
	 x0 -= wcsref[img][0];
	 y0 -= wcsref[img][1];
	 // apply SIP corrections
	 *x=x0
	    +x0*x0*wcs[img][2]+x0*y0*wcs[img][3]+y0*y0*wcs[img][4]
	    +x0*x0*x0*wcs[img][5]+x0*x0*y0*wcs[img][6]+x0*y0*y0*wcs[img][7]+y0*y0*y0*wcs[img][8]
	    +x0*x0*x0*x0*wcs[img][9]+x0*x0*x0*y0*wcs[img][10]+x0*x0*y0*y0*wcs[img][11]+x0*y0*y0*y0*wcs[img][12]+y0*y0*y0*y0*wcs[img][13]
	    +x0*x0*x0*x0*x0*wcs[img][14]+x0*x0*x0*x0*y0*wcs[img][15]+x0*x0*x0*y0*y0*wcs[img][16]+x0*x0*y0*y0*y0*wcs[img][17]+x0*y0*y0*y0*y0*wcs[img][18]+y0*y0*y0*y0*y0*wcs[img][19];
	 *y=y0
	    +x0*x0*wcs[img][22]+x0*y0*wcs[img][23]+y0*y0*wcs[img][24]
	    +x0*x0*x0*wcs[img][25]+x0*x0*y0*wcs[img][26]+x0*y0*y0*wcs[img][27]+y0*y0*y0*wcs[img][28]
	    +x0*x0*x0*x0*wcs[img][29]+x0*x0*x0*y0*wcs[img][30]+x0*x0*y0*y0*wcs[img][31]+x0*y0*y0*y0*wcs[img][32]+y0*y0*y0*y0*wcs[img][33]
	    +x0*x0*x0*x0*x0*wcs[img][34]+x0*x0*x0*x0*y0*wcs[img][35]+x0*x0*x0*y0*y0*wcs[img][36]+x0*x0*y0*y0*y0*wcs[img][37]+x0*y0*y0*y0*y0*wcs[img][38]+y0*y0*y0*y0*y0*wcs[img][39];
	 // apply 1st order corrections
	 x0=(*x)*wcs[img][0]+(*y)*wcs[img][1];
	 y0=(*x)*wcs[img][20]+(*y)*wcs[img][21];
      }
      else
#ifdef USEWFPC2
      if (hstmode[img].inst==WFPC2 && img<Nimg && DRIZZLE_BASE) {
	 x0+=hstoffset[img][0];
	 y0+=hstoffset[img][1];
	 WFPC2shift(img,&x0,&y0);
      }
      else
#endif
#ifdef USEACS
      if (hstmode[img].inst==ACS && img<Nimg && DRIZZLE_BASE) {
	 x0+=hstoffset[img][0];
	 y0+=hstoffset[img][1];
	 ACSshift(img,&x0,&y0);
      }
      else
#endif
#ifdef USEWFC3
      if (hstmode[img].inst==WFC3 && img<Nimg && DRIZZLE_BASE) {
	 x0+=hstoffset[img][0];
	 y0+=hstoffset[img][1];
	 WFC3shift(img,&x0,&y0);
      }
      else
#endif
#ifdef USEROMAN
      if (hstmode[img].inst==ROMAN && img<Nimg && DRIZZLE_BASE) {
	 x0+=hstoffset[img][0];
	 y0+=hstoffset[img][1];
	 ROMANshift(img,&x0,&y0);
      }
      else
#endif
#ifdef USENIRCAM
      if (hstmode[img].inst==NIRCAM && img<Nimg && DRIZZLE_BASE) {
	 x0+=hstoffset[img][0];
	 y0+=hstoffset[img][1];
	 NIRCAMshift(img,&x0,&y0);
      }
      else
#endif
#ifdef USENIRISS
      if (hstmode[img].inst==NIRISS && img<Nimg && DRIZZLE_BASE) {
	 x0+=hstoffset[img][0];
	 y0+=hstoffset[img][1];
	 NIRISSshift(img,&x0,&y0);
      }
      else
#endif
#ifdef USEMIRI
      if (hstmode[img].inst==MIRI && img<Nimg && DRIZZLE_BASE) {
	 x0+=hstoffset[img][0];
	 y0+=hstoffset[img][1];
	 MIRIshift(img,&x0,&y0);
      }
      else
#endif
      {
	 x0-=X*0.5;
	 y0-=Y*0.5;
      }
      if (dir==-1 || dir==-2) {
	 *x=x0*ref2img[img][21]+y0*ref2img[img][22]+x0*x0*ref2img[img][23]+x0*y0*ref2img[img][24]+y0*y0*ref2img[img][25]+x0*x0*x0*ref2img[img][26]+x0*x0*y0*ref2img[img][27]+x0*y0*y0*ref2img[img][28]+y0*y0*y0*ref2img[img][29];
	 *y=x0*ref2img[img][31]+y0*ref2img[img][32]+x0*x0*ref2img[img][33]+x0*y0*ref2img[img][34]+y0*y0*ref2img[img][35]+x0*x0*x0*ref2img[img][36]+x0*x0*y0*ref2img[img][37]+x0*y0*y0*ref2img[img][38]+y0*y0*y0*ref2img[img][39];
      }
      else {
	 *x=x0;
	 *y=y0;
      }
      if (dir==-1) {
	 x0 = *x; y0 = *y;
	 if (Align<4) {
	    x0 -= dpos[img][0];
	    y0 -= dpos[img][1];
	    r2=dpos[img][2]*(0.75+0.25*sqrt(1.+8*d2[img]*(x0*x0+y0*y0)/dpos[img][2]/dpos[img][2]));
	    *x=(c[img]*x0+s[img]*y0)/r2;
	    *y=(c[img]*y0-s[img]*x0)/r2;
	 }
	 else {
	    *x=dpos[img][20]+x0*dpos[img][21]+y0*dpos[img][22]+x0*x0*dpos[img][23]+x0*y0*dpos[img][24]+y0*y0*dpos[img][25]+x0*x0*x0*dpos[img][26]+x0*x0*y0*dpos[img][27]+x0*y0*y0*dpos[img][28]+y0*y0*y0*dpos[img][29];
	    *y=dpos[img][30]+x0*dpos[img][31]+y0*dpos[img][32]+x0*x0*dpos[img][33]+x0*y0*dpos[img][34]+y0*y0*dpos[img][35]+x0*x0*x0*dpos[img][36]+x0*x0*y0*dpos[img][37]+x0*y0*y0*dpos[img][38]+y0*y0*y0*dpos[img][39];
	 }
      }
      if (UseWCS==2) {
	 (*x)+=wcsref[img][2];
	 (*y)+=wcsref[img][3];
      }
      else {
	 (*x)+=X*0.5;
	 (*y)+=Y*0.5;
      }
   }

   // reverse transformation (reference to image)
   else {
      if (UseWCS==2) {
	 if (img<Nimg && wcsminmax_set && (x0<wcsminmax[img][0] || x0>wcsminmax[img][1] || y0<wcsminmax[img][2] || y0>wcsminmax[img][3])) {
	    *y = -1.e3;
	    *x = -1.e3;
	    return;
	 }
	 x0-=wcsref[img][2];
	 y0-=wcsref[img][3];
      }
      else {
	 x0-=X*0.5;
	 y0-=Y*0.5;
      }
      if (dir==1) {
	 if (Align<4) {
	    r2=dpos[img][2]*(1+d2[img]*(x0*x0+y0*y0));
	    *x=dpos[img][0]+(c[img]*x0-s[img]*y0)*r2;
	    *y=dpos[img][1]+(s[img]*x0+c[img]*y0)*r2;
	 }
	 else {
	    *x=dpos[img][0]+x0*dpos[img][1]+y0*dpos[img][2]+x0*x0*dpos[img][3]+x0*y0*dpos[img][4]+y0*y0*dpos[img][5]+x0*x0*x0*dpos[img][6]+x0*x0*y0*dpos[img][7]+x0*y0*y0*dpos[img][8]+y0*y0*y0*dpos[img][9];
	    *y=dpos[img][10]+x0*dpos[img][11]+y0*dpos[img][12]+x0*x0*dpos[img][13]+x0*y0*dpos[img][14]+y0*y0*dpos[img][15]+x0*x0*x0*dpos[img][16]+x0*x0*y0*dpos[img][17]+x0*y0*y0*dpos[img][18]+y0*y0*y0*dpos[img][19];
	 }
	 x0 = *x; y0 = *y;
      }
      if (dir==1 || dir==2) {
	 *x=x0*ref2img[img][1]+y0*ref2img[img][2]+x0*x0*ref2img[img][3]+x0*y0*ref2img[img][4]+y0*y0*ref2img[img][5]+x0*x0*x0*ref2img[img][6]+x0*x0*y0*ref2img[img][7]+x0*y0*y0*ref2img[img][8]+y0*y0*y0*ref2img[img][9];
	 *y=x0*ref2img[img][11]+y0*ref2img[img][12]+x0*x0*ref2img[img][13]+x0*y0*ref2img[img][14]+y0*y0*ref2img[img][15]+x0*x0*x0*ref2img[img][16]+x0*x0*y0*ref2img[img][17]+x0*y0*y0*ref2img[img][18]+y0*y0*y0*ref2img[img][19];
      }
      else {
	 *x=x0;
	 *y=y0;
      }
      if (UseWCS==2) {
	 // apply SIP correction
	 x0=(*x)
	    +(*x)*(*x)*wcs[img][42]+(*x)*(*y)*wcs[img][43]+(*y)*(*y)*wcs[img][44]
	    +(*x)*(*x)*(*x)*wcs[img][45]+(*x)*(*x)*(*y)*wcs[img][46]+(*x)*(*y)*(*y)*wcs[img][47]+(*y)*(*y)*(*y)*wcs[img][48]
	    +(*x)*(*x)*(*x)*(*x)*wcs[img][49]+(*x)*(*x)*(*x)*(*y)*wcs[img][50]+(*x)*(*x)*(*y)*(*y)*wcs[img][51]+(*x)*(*y)*(*y)*(*y)*wcs[img][52]+(*y)*(*y)*(*y)*(*y)*wcs[img][53]
	    +(*x)*(*x)*(*x)*(*x)*(*x)*wcs[img][54]+(*x)*(*x)*(*x)*(*x)*(*y)*wcs[img][55]+(*x)*(*x)*(*x)*(*y)*(*y)*wcs[img][56]+(*x)*(*x)*(*y)*(*y)*(*y)*wcs[img][57]+(*x)*(*y)*(*y)*(*y)*(*y)*wcs[img][58]+(*y)*(*y)*(*y)*(*y)*(*y)*wcs[img][59];
	 y0=(*y)
	    +(*x)*(*x)*wcs[img][62]+(*x)*(*y)*wcs[img][63]+(*y)*(*y)*wcs[img][64]
	    +(*x)*(*x)*(*x)*wcs[img][65]+(*x)*(*x)*(*y)*wcs[img][66]+(*x)*(*y)*(*y)*wcs[img][67]+(*y)*(*y)*(*y)*wcs[img][68]
	    +(*x)*(*x)*(*x)*(*x)*wcs[img][69]+(*x)*(*x)*(*x)*(*y)*wcs[img][70]+(*x)*(*x)*(*y)*(*y)*wcs[img][71]+(*x)*(*y)*(*y)*(*y)*wcs[img][72]+(*y)*(*y)*(*y)*(*y)*wcs[img][73]
	    +(*x)*(*x)*(*x)*(*x)*(*x)*wcs[img][74]+(*x)*(*x)*(*x)*(*x)*(*y)*wcs[img][75]+(*x)*(*x)*(*x)*(*y)*(*y)*wcs[img][76]+(*x)*(*x)*(*y)*(*y)*(*y)*wcs[img][77]+(*x)*(*y)*(*y)*(*y)*(*y)*wcs[img][78]+(*y)*(*y)*(*y)*(*y)*(*y)*wcs[img][79];
	 // apply 1st order correction and add reference pixel
	 *x=wcsref[img][0]+x0*wcs[img][40]+y0*wcs[img][41];
	 *y=wcsref[img][1]+x0*wcs[img][60]+y0*wcs[img][61];
      }
      else
#ifdef USEWFPC2
      if (hstmode[img].inst==WFPC2 && img<Nimg && DRIZZLE_BASE) {
	 WFPC2unshift(img,x,y);
	 (*x)-=hstoffset[img][0];
	 (*y)-=hstoffset[img][1];
      }
      else
#endif
#ifdef USEACS
      if (hstmode[img].inst==ACS && img<Nimg && DRIZZLE_BASE) {
	 ACSunshift(img,x,y);
	 (*x)-=hstoffset[img][0];
	 (*y)-=hstoffset[img][1];
      }
      else
#endif
#ifdef USEWFC3
      if (hstmode[img].inst==WFC3 && img<Nimg && DRIZZLE_BASE) {
	 WFC3unshift(img,x,y);
	 (*x)-=hstoffset[img][0];
	 (*y)-=hstoffset[img][1];
      }
      else
#endif
#ifdef USEROMAN
      if (hstmode[img].inst==ROMAN && img<Nimg && DRIZZLE_BASE) {
	 ROMANunshift(img,x,y);
	 (*x)-=hstoffset[img][0];
	 (*y)-=hstoffset[img][1];
      }
      else
#endif
#ifdef USENIRCAM
      if (hstmode[img].inst==NIRCAM && img<Nimg && DRIZZLE_BASE) {
	 NIRCAMunshift(img,x,y);
	 (*x)-=hstoffset[img][0];
	 (*y)-=hstoffset[img][1];
      }
      else
#endif
#ifdef USENIRISS
      if (hstmode[img].inst==NIRISS && img<Nimg && DRIZZLE_BASE) {
	 NIRISSunshift(img,x,y);
	 (*x)-=hstoffset[img][0];
	 (*y)-=hstoffset[img][1];
      }
      else
#endif
#ifdef USEMIRI
      if (hstmode[img].inst==MIRI && img<Nimg && DRIZZLE_BASE) {
	 MIRIunshift(img,x,y);
	 (*x)-=hstoffset[img][0];
	 (*y)-=hstoffset[img][1];
      }
      else
#endif
      {
	 (*x)+=X*0.5;
	 (*y)+=Y*0.5;
      }
   }
   return;
}

void setwcsminmax(void)
{
   int img,x0,y0;
   double x0f,y0f,x,y;
   for (img=0;img<Nimg;img++) {
      for (x0=0;x0<=16;x0++) {
	 // x0f/y0f code pads region so stars within RPSF of the image boundary will can have positions transformed
	 if (x0==0) x0f = -RPSF[img];
	 else if (x0==16) x0f = dataim[img].X + RPSF[img];
	 else x0f = x0*dataim[img].X/16.0;
	 for (y0=0;y0<=16;y0++) {
	    if (y0==0) y0f = -RPSF[img];
	    else if (y0==16) y0f = dataim[img].Y + RPSF[img];
	    else y0f = y0*dataim[img].Y/16.0;
	    //shift(img,x0*dataim[img].X/16.0,y0*dataim[img].Y/16.0,&x,&y,-1);
	    shift(img,x0f,y0f,&x,&y,-1);
	    if (x0==0 && y0==0) {
	       wcsminmax[img][0] = wcsminmax[img][1] = x;
	       wcsminmax[img][2] = wcsminmax[img][3] = y;
	    }
	    else {
	       if (x<wcsminmax[img][0]) wcsminmax[img][0] = x;
	       else if (x>wcsminmax[img][1]) wcsminmax[img][1] = x;
	       if (y<wcsminmax[img][2]) wcsminmax[img][2] = y;
	       else if (y>wcsminmax[img][3]) wcsminmax[img][3] = y;
	    }
	 }
      }
      //printf("%g-%g, %g-%g\n",wcsminmax[img][0],wcsminmax[img][1],wcsminmax[img][2],wcsminmax[img][3]);
      wcsminmax[img][0]-=5;
      wcsminmax[img][1]+=5;
      wcsminmax[img][2]-=5;
      wcsminmax[img][3]+=5;
      //printf("%g-%g, %g-%g\n",wcsminmax[img][0],wcsminmax[img][1],wcsminmax[img][2],wcsminmax[img][3]);
   }
   wcsminmax_set = 1;
}

void fixRA(int img,double ra0,double dec0,double *ra,double *dec,int dir)
{
   // from http://lambda.gsfc.nasa.gov/product/iras/coordproj.cfm

   double dra,ddec;
   dra = ra0-wcsref[img][2];
   ddec = dec0-wcsref[img][3];
   // forward correction: RA-TAN to RA
   if (dir<0) {
#if 1
      double X,Y,D,B,XX,YY;
      X = -dra*M_PI/180.0;
      Y = -ddec*M_PI/180.0;
      D = atan(sqrt(X*X+Y*Y)); // angle with dY=distance in radians, dX=1
      B = atan2(-X,Y); // angle from south
      XX = sin(wcsref[img][3]*M_PI/180.0)*sin(D)*cos(B) + cos(wcsref[img][3]*M_PI/180.0)*cos(D);
      YY = sin(D) * sin(B);
      // *ra ~= ra_c + (ra0-ra_c) / (cos(dec_c)-(dec0-dec_c)*sin(dec_c))
      *ra = wcsref[img][2] + atan2(YY,XX)*180/M_PI;
      // *dec ~= asin( sin(dec_c) + (dec0-dec_c)*cos(dec_c) )
      *dec = asin( sin(wcsref[img][3]*M_PI/180.0)*cos(D) - cos(wcsref[img][3]*M_PI/180.0)*sin(D)*cos(B) )*180/M_PI;
#else
      double RASF,DecSF;
      RASF = cos(dec0*M_PI/180.0);
      DecSF = sin(wcsref[img][3]*M_PI/180.0);
      *dec = wcsref[img][3] + ddec - (1-cos(dra/RASF*M_PI/180.0))*180/M_PI*RASF*DecSF;
      RASF = cos((*dec)*M_PI/180.0); // iterate
      *ra = wcsref[img][2] + dra/RASF;
      *dec = wcsref[img][3] + ddec - (1-cos(dra/RASF*M_PI/180.0))*180/M_PI*RASF*DecSF;
#endif
   }
   // reverse correction: RA to RA-TAN
   else {
#if 1
      double A,F;
      A = cos(dec0*M_PI/180.0) * cos(dra*M_PI/180.0);
      // Note F = (180.0/M_PI) / ( cos(dec-dec_c) - cos(dec_c)*cos(dec)*(1-cos(ra-ra_c)) ) -- basically equal to 180/M_PI for any DOLPHOT application
      F = (180.0/M_PI) / ( sin(wcsref[img][3]*M_PI/180.0) * sin(dec0*M_PI/180.0) + A * cos(wcsref[img][3]*M_PI/180.0) );
      *ra = wcsref[img][2] + F * cos(dec0*M_PI/180.0) * sin((ra0-wcsref[img][2])*M_PI/180.0);
      // Note *dec = dec_c + F * [ sin(dec-dec_c) + sin(dec_c)*cos(dec)*(1-cos(ra-ra_c)) ] -- first term is linear mapping, second is curvature of dec lines
      *dec = wcsref[img][3] + F * ( cos(wcsref[img][3]*M_PI/180.0) * sin(dec0*M_PI/180.0) - A * sin(wcsref[img][3]*M_PI/180.0) );
#else
      double RASF,DecSF;
      RASF = cos(dec0*M_PI/180.0);
      DecSF = sin(wcsref[img][3]*M_PI/180.0);
      *ra = wcsref[img][2] + dra*RASF;
      *dec = wcsref[img][3] + ddec + (1-cos(dra*M_PI/180.0))*180/M_PI*RASF*DecSF;
#endif
   }
}

// from Numerical Recipes
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(void)
{
   static long idum=0;
   int j;
   long k;
   static long idum2=123456789;
   static long iy=0;
   static long iv[NTAB];
   float temp;

   if (idum <= 0 || MCSEED > 0) {
      if (MCSEED > 0) idum = MCSEED;
      else idum = time(NULL);
      MCSEED = 0;
      idum2=idum;
      for (j=NTAB+7;j>=0;j--) {
	 k=idum/IQ1;
	 idum=IA1*(idum-k*IQ1)-k*IR1;
	 if (idum < 0) idum += IM1;
	 if (j < NTAB) iv[j] = idum;
      }
      iy=iv[0];
   }
   k=idum/IQ1;
   idum=IA1*(idum-k*IQ1)-k*IR1;
   if (idum < 0) idum += IM1;
   k=idum2/IQ2;
   idum2=IA2*(idum2-k*IQ2)-k*IR2;
   if (idum2 < 0) idum2 += IM2;
   j=iy/NDIV;
   iy=iv[j]-idum2;
   iv[j] = idum;
   if (iy < 1) iy += IMM1;
   if ((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/*
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double ran1(void) {
   int j;
   long k;
   static long idum=0;
   static long iy=0;
   static long iv[NTAB];
   double temp;

   if (idum <= 0 || !iy) {
      idum=time(NULL);
      if (idum < 1) idum=1;
      for (j=NTAB+7;j>=0;j--) {
	 k=idum/IQ;
	 idum=IA*(idum-k*IQ)-IR*k;
	 if (idum < 0) idum += IM;
	 if (j < NTAB) iv[j] = idum;
      }
      iy=iv[0];
   }
   k=idum/IQ;
   idum=IA*(idum-k*IQ)-IR*k;
   if (idum < 0) idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j] = idum;
   if ((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}	
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
*/

double gauss(void) {
   static int iset=0;
   static double gset;
   double fac,rsq,v1,v2;

   if  (iset == 0) {
      do {
	 v1=2.0*ran2()-1.0;
	 v2=2.0*ran2()-1.0;
	 rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return v2*fac;
   } else {
      iset=0;
      return gset;
   }
}

int poiss2(double m) {
   int n,mmin,mmax;
   double logm,pmax;

   if (m<=0.) {
      n=-(int)(-m);
      if (ran2()<n-m) n--;
      return n;
   }
   // set 5-sigma bounds
   mmin=(int)(m-5*sqrt(m)-1);
   if (mmin<0) mmin=0;
   mmax=(int)(m+5*sqrt(m)+2)-mmin+1;
   // P(n) = e^-m m^n / n!
   logm=log(m);
   n = (int)m;
   pmax = n*logm-lgamma(n+1);
   do {
      n=mmin+(int)(mmax*ran2());
   } while (ran2()>=exp(n*logm-lgamma(n+1)-pmax));
   return n;
}

double poiss3(double m) {
   if (m<=0.) {
      int n=-(int)(-m);
      if (ran2()<n-m) n--;
      return n;
   }
   return m + sqrt(m)*gauss();
}

int poiss1(double m) {
   int n=0;
   if (m<=0.) {
      n=-(int)(-m);
      if (ran2()<n-m) n--;
      return n;
   }
   if (m>=100) return poiss2(m);
   double p1 = exp(-m);
   double p = p1;
   double P = ran2();
   while (p<P) {
      n++;
      p1 *= m/n;
      p += p1;
   }
   return n;
}

int poiss(double m) {
   if (m<85) return poiss1(m);
   if (m<1.e6) poiss2(m);
   return poiss3(m);
}
