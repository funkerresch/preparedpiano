#include "ext.h"
#include "ext_obex.h"	// required for new style Max object
#include "ext_common.h"
#include "z_dsp.h"
#include "math.h"

#include "t_std.h"
#include "t_pianohammer.h"
#include "t_rattle.h"
#include "t_rubber.h"
#include "t_trap.h"

void *ppiano_class;
																									
/*																											    */
/*-------------------------------------------struct ppiano------------------------------------------------------*/
/*			  																									*/

typedef struct _ppiano
{
    t_pxobject x_obj;
    t_pianoHammer **hammer;
    t_rattle **rattle;
    t_rubber **rubber;
    t_trap **trap;
    
    memError err;
    
    double frequency;

    double **string;
    double **lastString;
    double **nextString;

    double dt;
    double dt_2;
    double dx;
    double dx_2;
    double dx_4;
    double f0;	//center freq
    double *f;	// freq strings
    double *c;	// c=T/p => wave speed = 2*f
    double *c_2;
    double *k;	// stiffness
    double *k_2;
    double *o;	// decay   [sigma]
    double *o_mul_dt_over_2;
    double *b;   // freq_dependent_loss
    double *alpha0;
    double *alpha1;
    double *alpha2;
    double *beta0;
    double *beta1;
    double *T30;   //  decay time (s)
    double *dxminAll;          
    double dxmin;
    double D;    // Detune in cents...
    double *detuneSpread;
    double scanFrequency;
    double scanFrequencySpread;
    double *listeningPointNorm;
    double scaleBack;
    double lastOut[10];
    float floatDummy;
    
    int *listeningPointIndex;
 
    int max_N;		// max Table length
    int N;			// Table Length
    int NS;  		//Number of Strings...
    int NH;		//Number of Hammers
    int numberOfPreparations;		//Number of rattles
    int NOutlets;
    int noteOff;
    
    char pitchAsMidi;
    char trap1;
    char trap2;
    char rattle1;
    char rattle2;
    char rubber1;
    char rubber2;
    
    char freqChanged;
    
    float trap1position;
    float trap2position;
    float rubber1position;
    float rubber2position;
    float rattle1position;
    float rattle2position;
    float rattle1freq;
    float rattle1massDensityRatio;
    float rattle1length;
    float rattle2freq;
    float rattle2massDensityRatio;
    float rattle2length;
    
    float attrDecay;
    float attrStiffness;
    float attrLoss;
    
} ppiano;

/*																												*/
/* ----------------------------------------- utility functions  ------------------------------------------------*/
/*																												*/


void setStringParametersHelp(ppiano *x, int count, int *stringNumber, double *decay, double *stiffness, double *freqDepLoss);
void ppiano_dsp64(ppiano *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
					
void calcCoefficients(ppiano *x, int stringNumber)
{
	int i = stringNumber;
	
    x->alpha0[i] = (2 -
                   ((2*x->c_2[i]*x->dt_2)/x->dx_2) -
                   ((6*x->k_2[i]*x->dt_2)/x->dx_4) -
                   ((2*x->b[i]*x->dt)/x->dx_2)) /
                   (1 + x->o_mul_dt_over_2[i]);
    
    x->alpha1[i] = (((x->c_2[i]*x->dt_2)/x->dx_2) +
                   ((4*x->k_2[i]*x->dt_2)/x->dx_4) +
                   ((x->b[i]*x->dt)/x->dx_2)) /
                   (1 + x->o_mul_dt_over_2[i]);
                        
    x->alpha2[i] = ((-x->k_2[i]*x->dt_2) /
                   (x->dx_4*(1 + x->o_mul_dt_over_2[i])));
                
    x->beta0[i] = (-1+2*x->b[i]*x->dt/x->dx_2+x->o_mul_dt_over_2[i]) /
                  (1 + x->o_mul_dt_over_2[i]);
               
    x->beta1[i] = (-x->b[i]*x->dt)/(x->dx_2 * (1 + x->o_mul_dt_over_2[i]));
}

void resetStrings(ppiano *x)
{
	int i;
	int m;
    
	for(i=0;i<x->NS;i++)
	{
 		for (m=0; m<=x->max_N; m++)
 		{
 			x->nextString[i][m] = 0;
 			x->string[i][m] = 0;
 			x->lastString[i][m] = 0;
 		}
    }
}

void resetAll(ppiano *x)
{
	int i;
	int m;
    
	for(i=0;i<x->NS;i++)
	{
 		for (m=0; m<=x->max_N; m++)
 		{
 			x->nextString[i][m] = 0;
 			x->string[i][m] = 0;
 			x->lastString[i][m] = 0;
 		}
    }
 
    t_pianoHammer_reset(x->hammer[0]);
    
    for(i=0;i<x->numberOfPreparations;i++)
    {
        t_rattle_reset(x->rattle[i]);	
        t_rubber_reset(x->rubber[i]);
        t_trap_reset(x->trap[i]);
    }
}

void setFrequency(ppiano *x, double frequency, double detune)
{
    int i;
    
	x->f0 = frequency;
	x->D = detune;
   	
	for(i=0;i<x->NS;i++)
	{	
		if(x->NS>1)
		{
			x->detuneSpread[i] = x->D*i/(x->NS-1) - x->D/2;
			x->f[i] = x->f0*pow(2, x->detuneSpread[i]/1200.);
		}
		else
			x->f[i] = x->f0;
		
		x->c[i] = 2*x->f[i];
		x->c_2[i] = pow(x->c[i], 2);
	}
	
	for(i=0;i<x->NS;i++)
		x->dxminAll[i] = sqrt(  ((x->c_2[i]*x->dt_2+2*x->b[i]*x->dt) + sqrt( pow(x->c_2[i]*x->dt_2+2*x->b[i]*x->dt,2)+16*x->k_2[i]*x->dt_2)) / 2.);
			
	x->dxmin = t_fmax(x->dxminAll, x->NS);
    
	x->N = floor(1/x->dxmin);
    
	x->dx = 1./(double)x->N;
	x->dx_2 = pow(x->dx,2);
	x->dx_4 = pow(x->dx,4);
	
	for(i=0;i<x->NOutlets;i++)
		x->listeningPointIndex[i] = 2+floor(x->listeningPointNorm[i]/x->dx);
	
	for(i=0;i<x->NS;i++)
		calcCoefficients(x, i);
    
    x->freqChanged = 1;
    
    t_pianoHammer_set_dxdt(x->hammer[0], x->dx, x->dt);
    
    for(i=0;i<2;i++)
    {
		t_rattle_set_dxdt(x->rattle[i], x->dx, x->dt);
		t_rubber_set_dxdt(x->rubber[i], x->dx, x->dt);
		t_trap_set_dxdt(x->trap[i], x->dx, x->dt);
    }
}

void setListeningPoint(ppiano *x, int outletNumber, double listeningPoint)
{
	if(listeningPoint >=0 && listeningPoint <=1 && outletNumber >= 0 && outletNumber < x->NOutlets)
	{
		x->listeningPointNorm[outletNumber] = listeningPoint;
		x->listeningPointIndex[outletNumber] = 2+floor(x->listeningPointNorm[outletNumber]/x->dx);
	}
}

void setHammerPosition(ppiano *x,  double hammerPosition)
{
	if(hammerPosition >=0 && hammerPosition <= 1)
	{
		x->hammer[0]->normHammerIndex = hammerPosition;
		t_pianoHammer_calculateHammerIndex(x->hammer[0]);
	}
	else
		post("error: hammerPosition must be >=0 and <= 1");
}

//-------------------------------------------------------Hammer-------------------------------------------------------------

void hammerStrike(ppiano *x, t_symbol *message, short argc, t_atom *argv)
{
    int hammerNumber = 1;
    double hammerVel = 80;
    
    if(argc >= 1)
    {
        if(atom_gettype(argv) == A_LONG)
            hammerVel = atom_getlong(argv);
        else if(atom_gettype(argv) == A_FLOAT)
            hammerVel = atom_getfloat(argv);
        else
            return;
    }
    
    if(hammerNumber <= x->NH)
    {
        int i;
        i = hammerNumber-1;
     
        x->hammer[i]->hammerVelocity = hammerVel;
        x->hammer[i]->lastHammer = x->hammer[i]->hammerInitPosition;
        x->hammer[i]->hammer = x->hammer[i]->hammerInitPosition + (x->dt) * x->hammer[i]->hammerVelocity;
        x->hammer[i]->nextHammer = 0;
        x->hammer[i]->hammerContact = 0;
        x->hammer[i]->hammerOn = 1;
    }
}

void hammerStrike1(ppiano *x, double velocity)
{
    double hammerVel = velocity;
        
    x->hammer[0]->hammerVelocity = hammerVel;
    x->hammer[0]->lastHammer = x->hammer[0]->hammerInitPosition;
    x->hammer[0]->hammer = x->hammer[0]->hammerInitPosition + (x->dt) * x->hammer[0]->hammerVelocity;
    x->hammer[0]->nextHammer = 0;
    x->hammer[0]->hammerContact = 0;
    x->hammer[0]->hammerOn = 1;
}

void plug1(ppiano *x, t_symbol *message, short argc, t_atom *argv)
{
    t_int32 i = 0;
    t_int32 plugIndex;
    
    double pitch = 72;
    double velocity = 80;
    double duration = 1;
    double plugPosition = 0.5;
    double detune = 0;
    
    if(argc < 2)
        return;
    
    x->noteOff = 0;
    
    if(atom_gettype(argv+i) == A_LONG)
        pitch = atom_getlong(argv+i);
    else if(atom_gettype(argv+i) == A_FLOAT)
        pitch = atom_getfloat(argv+i);
    else return;
    
    i++;
    
    if(atom_gettype(argv+i) == A_LONG)
        velocity = atom_getlong(argv+i);
    else if(atom_gettype(argv+i) == A_FLOAT)
        velocity = atom_getfloat(argv+i);
    else return;
    
    i++;
    
    if(i < argc)
    {
        if(atom_gettype(argv+i) == A_LONG)
            duration = atom_getlong(argv+i);
        else if(atom_gettype(argv+i) == A_FLOAT)
            duration = atom_getfloat(argv+i);
        else return;
        
    i++;
    }
    
    
    if(i < argc)
    {
        if(atom_gettype(argv+i) == A_LONG)
            plugPosition = atom_getlong(argv+i);
        else if(atom_gettype(argv+i) == A_FLOAT)
            plugPosition = atom_getfloat(argv+i);
        else return;
        
        i++;
    }
    
    if(i < argc)
    {
        if(atom_gettype(argv+i) == A_LONG)
            detune = atom_getlong(argv+i);
        else if(atom_gettype(argv+i) == A_FLOAT)
            detune = atom_getfloat(argv+i);
        else return;
        
        i++;
    }
    
    if(object_attr_get(x, gensym("pitchasmidinote")))
        pitch = midi2freq(pitch, 440);
    
    setFrequency(x, pitch, detune);
    plugIndex= 2+floor(plugPosition/x->dx);
    
    for(i = 0; i<x->NS; i++)
        x->string[i][plugIndex] = velocity / 500000.;
}

void noteOff(ppiano *x)
{
    x->noteOff = 1;
}

void noteon(ppiano *x, t_symbol *message, short argc, t_atom *argv)
{
    t_int32 i = 0;
   
    double pitch = 72;
    double velocity = 80;
    double duration = 1;
    double hammerPosition = 0.5;
    double detune = 0;
    x->noteOff = 0;
    
    if(argc < 2)
        return;
    
    if(atom_gettype(argv+i) == A_LONG)
        pitch = atom_getlong(argv+i);
    else if(atom_gettype(argv+i) == A_FLOAT)
        pitch = atom_getfloat(argv+i);
    else return;
    
    i++;
    
    if(atom_gettype(argv+i) == A_LONG)
        velocity = atom_getlong(argv+i);
    else if(atom_gettype(argv+i) == A_FLOAT)
        velocity = atom_getfloat(argv+i);
    else return;
    
    i++;
    
    if(i < argc)
    {
        if(atom_gettype(argv+i) == A_LONG)
            duration = atom_getlong(argv+i);
        else if(atom_gettype(argv+i) == A_FLOAT)
            duration = atom_getfloat(argv+i);
        else return;
        
         i++;
    }
    
   
    
    if(i < argc)
    {
        if(atom_gettype(argv+i) == A_LONG)
            hammerPosition = atom_getlong(argv+i);
        else if(atom_gettype(argv+i) == A_FLOAT)
            hammerPosition = atom_getfloat(argv+i);
        else return;
        
        i++;
    }
    
    if(i < argc)
    {
        if(atom_gettype(argv+i) == A_LONG)
            detune = atom_getlong(argv+i);
        else if(atom_gettype(argv+i) == A_FLOAT)
            detune = atom_getfloat(argv+i);
        else return;
        
        i++;
    }
    
    if(object_attr_get(x, gensym("pitchasmidinote")))
        pitch = midi2freq(pitch, 440);
    
    setHammerPosition(x, hammerPosition);
    setFrequency(x, pitch, detune);
    hammerStrike1(x, velocity);
}

//------------------------------------------------------Rattle----------------------------------------------------------------

void setRattlePosition(ppiano *x, int rattleNumber, double rattlePosition)
{
	int i;
	i = rattleNumber-1;
	if(rattlePosition >=0 && rattlePosition <= 1)
    {
		t_rattle_calculateRattleIndex(x->rattle[i], rattlePosition);
    }
	else
		post("error: rattlePosition must be >=0 and <= 1");
}

void setTrapFrequency(ppiano *x, int trapNumber, double trapFrequency)
{
    int i;
    i = trapNumber-1;
    t_trap_setFundamentalFrequency(x->trap[i], trapFrequency);
}

void setTrapMassDensityRation(ppiano *x, int trapNumber, double massDensityRatio)
{
    int i;
    i = trapNumber-1;
    t_trap_setMassDensityRatio(x->trap[i], massDensityRatio);
}

void setRattleFrequency(ppiano *x, int rattleNumber, double rattleFrequency)
{
    int i;
    i = rattleNumber-1;
    t_rattle_setFundamentalFrequency(x->rattle[i], rattleFrequency);
}

void setRattleMassDensityRatio(ppiano *x, int rattleNumber, double massDensityRation)
{
    int i;
    i = rattleNumber-1;
    t_rattle_setMassDensityRatio(x->rattle[i], massDensityRation);
}

void setRattleLength(ppiano *x, int rattleNumber, double rattleLength)
{
    int i;
    i = rattleNumber-1;
    t_rattle_setRattleLength(x->rattle[i], rattleLength);
}

void rattleOn(ppiano *x, int rattleNumber, int on)
{
	if(rattleNumber <= x->numberOfPreparations)
    {
   		int i;
   		i = rattleNumber-1;	
   		x->rattle[i]->rattleOn = on;
    }
}



//------------------------------------------------------- Rubber ----------------------------------------------------------------

void setRubberPosition(ppiano *x, int rubberNumber, double rubberPosition)
{
	int i;
	i = rubberNumber-1;
	if(rubberPosition >=0 && rubberPosition <= 1)
		t_rubber_calculateRubberIndex(x->rubber[i], rubberPosition);
	else
		post("error: rubberPosition must be >=0 and <= 1");
}

void rubberOn(ppiano *x, int rubberNumber, int on)
{
	if(rubberNumber <= x->numberOfPreparations)
    {
   		int i;
   		i = rubberNumber-1;	
   		x->rubber[i]->rubberOn = on;
    }
}

void setRubberParameters(ppiano *x, t_symbol *message, short argc, t_atom *argv)     //rubberFundamentalFreq, massDensityRatio, rubberLoss
{
	int i;
	int j;
	int s;
	int t;
	double k;
	int *rubberNumber;
	double *fundamentalFreq;
	double *massDensityRatio;
	double *rubberLoss;
	
	if(argc < 2)
		post("error: at least 2 arguments must be specified, first argument is rubber-number, second decay-rate");
	else 
	{
		k = (double)(argc)/4.;
		s = ceil(k);
		rubberNumber = (int *)sysmem_newptrclear(sizeof(int)*s);
		fundamentalFreq = (double *)sysmem_newptrclear(sizeof(double)*s);
		massDensityRatio = (double *)sysmem_newptrclear(sizeof(double)*s);
		rubberLoss = (double *)sysmem_newptrclear(sizeof(double)*s);
		for(t=0;t<s;t++)
		{
			rubberNumber[t] = -1;
			fundamentalFreq[t] = -1;
			massDensityRatio[t] = -1;
			rubberLoss[t] = -1;
		}
	
		for(i=0;i<argc;i++)
		{
			j=i%4;
			t=i/4;
			if(j==0)
			{
				switch (argv[i].a_type) 
				{        
					case A_LONG:  
						if(argv[i].a_w.w_long>x->numberOfPreparations || argv[i].a_w.w_long<=0)
						{
							post("error: selected Rubber doesn't exist....");
							goto SETPARAMEND;
						}
						else
							rubberNumber[t] = argv[i].a_w.w_long;
						break;        
					case A_SYM:         
						post("error: rubber-number must be an integer");         
						break;        
					case A_FLOAT:          
						post("error: rubber-number must be an integer");         
						break;      
		       
				 }  
			}
			else if(j==1)
			{
				 switch (argv[i].a_type) 
				 {        
					case A_LONG:  
						fundamentalFreq[t] = argv[i].a_w.w_long;
						break;        
					case A_SYM:         
						post("error: fundamentalFreq must be an integer or double");         
						break;        
					case A_FLOAT:          
						fundamentalFreq[t] = argv[i].a_w.w_float;        
						break;      
		       
				 }  
			}
			else if(j==2)
			{
				 switch (argv[i].a_type) 
				 {        
					case A_LONG:  
						massDensityRatio[t] = argv[i].a_w.w_long;
						break;        
					case A_SYM:         
						post("error: massDensityRatio must be an integer or double");         
						break;        
					case A_FLOAT:          
						massDensityRatio[t] = argv[i].a_w.w_float;        
						break;      
		       
				 }  
			}
			else if(j==3)
			{
				 switch (argv[i].a_type) 
				 {        
					case A_LONG:  
						rubberLoss[t] = argv[i].a_w.w_long;
						break;        
					case A_SYM:         
						post("error: rubberLoss must be an integer or double");         
						break;        
					case A_FLOAT:          
						rubberLoss[t] = argv[i].a_w.w_float;        
						break;      
			     }  
			}
	  
		}

		for(i=0;i<s;i++)
		{
			j=rubberNumber[i]-1;
			t_rubber_setRubberParameters(x->rubber[j], fundamentalFreq[i], massDensityRatio[i], rubberLoss[i]);
		}
SETPARAMEND:
        sysmem_freeptr(rubberNumber);
        sysmem_freeptr(fundamentalFreq);
        sysmem_freeptr(massDensityRatio);
        sysmem_freeptr(rubberLoss);
  			// sorry, but in this case the label/goto-statement does make sense......:)
	}
}

//---------------------------------------------------------Trap----------------------------------------------------------------

void trapOn(ppiano *x, int trapNumber, int on)
{
	if(trapNumber <= x->numberOfPreparations)
	{
		int i;
		i = trapNumber-1;
		x->trap[i]->trapOn = on;	
	}
}

void setTrapPosition(ppiano *x, int trapNumber, double trapPosition)
{
	int i;
	i = trapNumber-1;
	if(trapPosition >=0 && trapPosition <= 1)
		t_trap_calculateTrapIndex(x->trap[i], trapPosition);
	else
		post("error: trapPosition must be >=0 and <= 1");
}

void setTrapParameters(ppiano *x, t_symbol *message, short argc, t_atom *argv)     //trapFundamentalFreq, massDensityRatio
{
	int i;
	int j;
	int s;
	int t;
	double k;
	int *trapNumber;
	double *fundamentalFreq;
	double *massDensityRatio;
	
	if(argc < 2)
		post("error: at least 2 arguments must be specified, first argument is trap-number, second fundamental-freq");
	else 
	{
		k = (double)(argc)/3.;
		s = ceil(k);
		trapNumber = (int *)sysmem_newptrclear(sizeof(int)*s);
		fundamentalFreq = (double *)sysmem_newptrclear(sizeof(double)*s);
		massDensityRatio = (double *)sysmem_newptrclear(sizeof(double)*s);
	
		for(t=0;t<s;t++)
		{
			trapNumber[t] = -1;
			fundamentalFreq[t] = -1;
			massDensityRatio[t] = -1;
		}
	
		for(i=0;i<argc;i++)
		{
			j=i%3;
			t=i/3;
			if(j==0)
			{
				switch (argv[i].a_type) 
				{        
					case A_LONG:  
						if(argv[i].a_w.w_long>x->numberOfPreparations || argv[i].a_w.w_long<=0)
						{
							post("error: selected Trap doesn't exist....");
							goto SETPARAMEND;
						}
						else
						{
							trapNumber[t] = argv[i].a_w.w_long;
					
						}
						
						break;        
					case A_SYM:         
						post("error: trap-number must be an integer");         
						break;        
					case A_FLOAT:          
						post("error: trap-number must be an integer");         
						break;      
		       
				 }  
			}
			else if(j==1)
			{
				 switch (argv[i].a_type) 
				 {        
					case A_LONG:  
						fundamentalFreq[t] = argv[i].a_w.w_long;
						break;        
					case A_SYM:         
						post("error: fundamentalFreq must be an integer or double");         
						break;        
					case A_FLOAT:          
						fundamentalFreq[t] = argv[i].a_w.w_float;        
						break;      
		       
				 }  
			}
			else if(j==2)
			{
				 switch (argv[i].a_type) 
				 {        
					case A_LONG:  
						massDensityRatio[t] = argv[i].a_w.w_long;
						break;        
					case A_SYM:         
						post("error: massDensityRatio must be an integer or double");         
						break;        
					case A_FLOAT:          
						massDensityRatio[t] = argv[i].a_w.w_float;        
						break;      
		       
				 }  
			}
		}

		for(i=0;i<s;i++)
		{
			j=trapNumber[i]-1;
			t_trap_setTrapParameters(x->trap[j], fundamentalFreq[i], massDensityRatio[i]);
		}
        
SETPARAMEND:
 	sysmem_freeptr(trapNumber);
 	sysmem_freeptr(fundamentalFreq);
	sysmem_freeptr(massDensityRatio);
  			
	}
}

void setTrapParameters1(ppiano *x, t_symbol *message, short argc, t_atom *argv)     //trapFundamentalFreq, massDensityRatio
{
    int trapNumber;
    double fundamentalFreq;
    double massDensityRatio;
    
    if(argc < 2)
        error("At least 2 arguments must be specified, first argument is trap-number, second fundamental-freq");
    else
    {
        if(atom_gettype(argv) == A_LONG)
            trapNumber = atom_getlong(argv);
        else
            return;
        
        if (atom_gettype(argv+1) == A_FLOAT)
        {
             fundamentalFreq = atom_getfloat(argv+1);
             setTrapFrequency(x, 1, fundamentalFreq);
        }
        else
            return;
        
        if(argc >= 3)
        {
            if (atom_gettype(argv+2) == A_FLOAT)
            {
                massDensityRatio = atom_getfloat(argv+2);
                setTrapMassDensityRation(x, trapNumber, massDensityRatio);
            }
            else
                return;
        }
    }
}

//--------------------------------------------------------String----------------------------------------------------------------

void setStringParameters(ppiano *x, t_symbol *message, short argc, t_atom *argv)
{
    int setAll = 0;
	int i;
	int j;
	int s;
	int t;
	double k;
	int *stringNumber;
	double *decay;
	double *stiffness;
	double *freqDepLoss;
	
	if(argc < 2)
		post("error: at least 2 arguments must be specified, first argument is string-number, second decay-rate");
	else 
	{
		k = (double)(argc)/4.;
		s = ceil(k);
		stringNumber = (int *)sysmem_newptrclear(sizeof(int)*s);
		decay = (double *)sysmem_newptrclear(sizeof(double)*s);
		stiffness = (double *)sysmem_newptrclear(sizeof(double)*s);
		freqDepLoss = (double *)sysmem_newptrclear(sizeof(double)*s);
		for(t=0;t<s;t++)
		{
			stringNumber[t] = -1;
			decay[t] = -1;
			stiffness[t] = -1;
			freqDepLoss[t] = -1;
		}
	
		for(i=0;i<argc;i++)
		{
			j=i%4;
			t=i/4;
			if(j==0)
			{
				switch (argv[i].a_type) 
				{        
					case A_LONG:  
						if(argv[i].a_w.w_long>x->NS || argv[i].a_w.w_long<=0)
						{
							post("error: selected String doesn't exist....");
							goto SETPARAMEND;
						}
						else
							stringNumber[t] = argv[i].a_w.w_long;
						break;        
					case A_SYM:         
						setAll = 1;        
						break;        
					case A_FLOAT:          
						post("error: string-number must be an integer");         
						break;      
				 }  
			}
			else if(j==1)
			{
				 switch (argv[i].a_type) 
				 {        
					case A_LONG:  
						decay[t] = argv[i].a_w.w_long;
						break;        
					case A_SYM:         
						post("error: decay_rate must be an integer or double");         
						break;        
					case A_FLOAT:          
						decay[t] = argv[i].a_w.w_float;        
						break;      
				 }  
			}
			else if(j==2)
			{
				 switch (argv[i].a_type) 
				 {        
					case A_LONG:  
						stiffness[t] = argv[i].a_w.w_long;
						break;        
					case A_SYM:         
						post("error: stiffness must be an integer or double");         
						break;        
					case A_FLOAT:          
						stiffness[t] = argv[i].a_w.w_float;        
						break;      
				 }  
			}
			else if(j==3)
			{
				 switch (argv[i].a_type) 
				 {        
					case A_LONG:  
						freqDepLoss[t] = argv[i].a_w.w_long;
						break;        
					case A_SYM:         
						post("error: frequencyDependentLoss must be an integer or double");         
						break;        
					case A_FLOAT:          
						freqDepLoss[t] = argv[i].a_w.w_float;        
						break;      
				 }  
			}
		}
	
        if(!setAll)
            setStringParametersHelp(x, s, stringNumber, decay, stiffness, freqDepLoss);
        
SETPARAMEND:
 	sysmem_freeptr(stringNumber);
 	sysmem_freeptr(decay);
	sysmem_freeptr(stiffness);
	sysmem_freeptr(freqDepLoss);
  			
	}
}

void setStringParameters1(ppiano *x, int j, double decay, double stiffness, double freqDepLoss) // j  is stringnumber
{
    int i=0;

    if(decay>0)
    {
        x->T30[j] = decay;
        x->o[j] = (2/x->dt)*(pow(10,(3*x->dt/x->T30[j]))-1);
        x->o_mul_dt_over_2[j] = x->o[j] * x->dt * 0.5;
    }
    if(stiffness>0)
    {
        x->k[j] = stiffness;
        x->k_2[j] = pow(x->k[j],2);
    }
    if(freqDepLoss>0)
        x->b[j] = freqDepLoss;
    
    for(i=0;i<x->NS;i++)
        x->dxminAll[i] = sqrt(  ((x->c_2[i]*x->dt_2+2*x->b[i]*x->dt) + sqrt( pow(x->c_2[i]*x->dt_2+2*x->b[i]*x->dt,2)+16*x->k_2[i]*x->dt_2)) / 2.);
    
    x->dxmin = t_fmax(x->dxminAll, x->NS);
    x->N = floor(1/x->dxmin);
    x->dx = 1./(double)x->N;
    x->dx_2 = pow(x->dx,2);
    x->dx_4 = pow(x->dx,4);
    
    for(i=0;i<x->NS;i++)
        calcCoefficients(x,i);
    
    x->hammer[0]->dx = x->dx;
    t_pianoHammer_calculateHammerIndex(x->hammer[0]);
    
    for(i=0;i<x->numberOfPreparations;i++)
    {
        x->rattle[i]->dx = x->dx;
        t_rattle_calculateRattleIndex(x->rattle[i], x->rattle[i]->normRattleIndex);
        x->rubber[i]->dx = x->dx;
        t_rubber_calculateRubberIndex(x->rubber[i], x->rubber[i]->normRubberIndex);
        x->trap[i]->dx = x->dx;
        t_trap_calculateTrapIndex(x->trap[i], x->trap[i]->normTrapIndex);
    }
    for(i=0;i<x->NOutlets;i++)
        x->listeningPointIndex[i] = 2+floor(x->listeningPointNorm[i]/x->dx);
}

void setStringDecay(ppiano *x, double decay)
{
    int i = 0;
    for(i = 0; i< x->NS; i++)
         setStringParameters1(x, i, decay, -1, -1);
}

void setStringStiffness(ppiano *x, double stiffness)
{
    int i = 0;
    for(i = 0; i< x->NS; i++)
        setStringParameters1(x, i, -1, stiffness, -1);
}

void setStringLoss(ppiano *x, double loss)
{
    int i = 0;
    for(i = 0; i< x->NS; i++)
        setStringParameters1(x, i, -1, -1, loss);
}

void setStringParametersHelp(ppiano *x, int count, int *stringNumber, double *decay, double *stiffness, double *freqDepLoss)
{
	int i;
	int j;
	for(i=0;i<count;i++)
	{
		j = stringNumber[i]-1;
		if(decay[j]>0)
		{
			x->T30[j] = decay[j];
			x->o[j] = (2/x->dt)*(pow(10,(3*x->dt/x->T30[j]))-1);
			x->o_mul_dt_over_2[j] = x->o[j] * x->dt * 0.5;
		}
		if(stiffness[j]>0)
		{
			x->k[j] = stiffness[j];
			x->k_2[j] = pow(x->k[j],2);
		}
		if(freqDepLoss[j]>0)
			x->b[j] = freqDepLoss[j];
	}
	for(i=0;i<x->NS;i++)
		x->dxminAll[i] = sqrt(  ((x->c_2[i]*x->dt_2+2*x->b[i]*x->dt) + sqrt( pow(x->c_2[i]*x->dt_2+2*x->b[i]*x->dt,2)+16*x->k_2[i]*x->dt_2)) / 2.);
	
	x->dxmin = t_fmax(x->dxminAll, x->NS);
	x->N = floor(1/x->dxmin);
	x->dx = 1./(double)x->N;
	x->dx_2 = pow(x->dx,2);
	x->dx_4 = pow(x->dx,4);
	
	for(i=0;i<x->NS;i++)
	{	
		calcCoefficients(x,i);
	}		

    x->hammer[0]->dx = x->dx;
    t_pianoHammer_calculateHammerIndex(x->hammer[0]);
	
	for(i=0;i<x->numberOfPreparations;i++)
	{
		x->rattle[i]->dx = x->dx;
		t_rattle_calculateRattleIndex(x->rattle[i], x->rattle[i]->normRattleIndex);
        x->rubber[i]->dx = x->dx;
		t_rubber_calculateRubberIndex(x->rubber[i], x->rubber[i]->normRubberIndex);
        x->trap[i]->dx = x->dx;
		t_trap_calculateTrapIndex(x->trap[i], x->trap[i]->normTrapIndex);
	}
	for(i=0;i<x->NOutlets;i++)
		x->listeningPointIndex[i] = 2+floor(x->listeningPointNorm[i]/x->dx);
}

void setNumberOfStrings(ppiano *x, int NS)
{
	if(NS>0)
	{
        int i;
        x->NS = NS;
        for(i=0;i<x->NS;i++)
        {
            if(x->NS>1)
            {
                x->detuneSpread[i] = x->D*i/(x->NS-1) - x->D/2;
                x->f[i] = x->f0*pow(2, x->detuneSpread[i]/1200.);
            }
            else
                x->f[i] = x->f0;
            
            x->T30[i] = 4;
            x->k[i] = 4;
            x->k_2[i] = pow(x->k[i],2);
            x->b[i] = 0.002;
            x->c[i] = 2*x->f[i];
            x->c_2[i] = pow(x->c[i], 2);
            x->o[i] = (2/x->dt)*(pow(10,(3*x->dt/x->T30[i]))-1);
            x->o_mul_dt_over_2[i] = x->o[i] * x->dt * 0.5;
            x->dxminAll[i] = sqrt(  ((x->c_2[i]*x->dt_2+2*x->b[i]*x->dt) + sqrt( pow(x->c_2[i]*x->dt_2+2*x->b[i]*x->dt,2)+16*x->k_2[i]*x->dt_2)) / 2.);
        }

        x->dxmin = t_fmax(x->dxminAll, x->NS);
        x->N = floor(1/x->dxmin);
        x->dx = 1./(double)x->N;
        x->dx_2 = pow(x->dx,2);
        x->dx_4 = pow(x->dx,4);
    
        for(i=0;i<x->NS;i++)
            calcCoefficients(x, i);
            
        for(i=0;i<x->NOutlets;i++)
            x->listeningPointIndex[i] = 2+floor(x->listeningPointNorm[i]/x->dx);

        for(i=0;i<x->NH;i++)
        {
            x->hammer[i]->NS = x->NS;
            t_pianoHammer_set_dxdt(x->hammer[i], x->dx, x->dt);
        }
        for(i=0;i<x->numberOfPreparations;i++)
        {
            x->rattle[i]->NS = x->NS;
            t_rattle_set_dxdt(x->rattle[i], x->dx, x->dt);
        
            x->rubber[i]->NS = x->NS;
            t_rubber_set_dxdt(x->rubber[i], x->dx, x->dt);

            x->trap[i]->NS = x->NS;
            t_trap_set_dxdt(x->trap[i], x->dx, x->dt);
        }
	}
	else
		post("error: NS must be greater than 0 and smaller than max_NS");
}

//																																							//
//----------------------------------------------------------------------------------------------------------------------------------------------------------//
//																																						    //
 
void *ppiano_new(t_symbol *s, short argc, t_atom *argv);
void ppiano_free(ppiano *x);

void *ppiano_new(t_symbol *s, short argc, t_atom *argv)
{
    ppiano *x = NULL;
    
    if((x = (ppiano *)object_alloc(ppiano_class)) )
    {
        dsp_setup((t_pxobject *)x, 1);
        int outlets = 2;
        int numberOfStrings = 1;
        int maxNS = 5;
        int i;
        x->noteOff = 1;
        
        if(argc>0)
        {
            for(i=0;i<argc;i++)
            {
                switch (argv[i].a_type) 
                {        
                    case A_LONG: 
                        if(i==0) 
                            outlets = argv[i].a_w.w_long;
                        else if(i==1)
                            numberOfStrings = argv[i].a_w.w_long;
                        else if(i==2)
                            maxNS = argv[i].a_w.w_long;
                        else	
                            ;
                        break;        
                    case A_SYM:         
                        post("error: string-number must be an integer");         
                        break;        
                    case A_FLOAT:          
                        post("error: string-number must be an integer");         
                        break;      
                }  
            }
        }
        
        if(outlets<=0)
            outlets = 1;
        if(numberOfStrings<=0)
            numberOfStrings = 1;
        if(maxNS<=0)
            maxNS = 3;
        if(maxNS<numberOfStrings)
            maxNS=numberOfStrings;
        
        x->freqChanged = 0;
        
        x->NOutlets = outlets;
        x->NS = numberOfStrings;

        for(i=1;i<=x->NOutlets;i++)
            outlet_new((t_pxobject *)x, "signal"); 
        
        if(numberOfStrings > 0)
            x->NS = numberOfStrings;
        
        x->max_N = 1000;
        x->NH = 1;
        x->numberOfPreparations = 2;
        x->f0 = 500;
        x->D = 1;
        x->scaleBack = 1;
            
        x->T30 = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->listeningPointNorm = (double *)sysmem_newptrclear(sizeof(double)*x->NOutlets);
        x->listeningPointIndex = (int *)sysmem_newptrclear(sizeof(int)*x->NOutlets);
        x->k = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->k_2 = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->b = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->o = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->o_mul_dt_over_2 = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->f = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->detuneSpread = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->c = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->c_2 = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->dxminAll = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->alpha0 = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->alpha1 = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->alpha2 = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->beta0 = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->beta1 = (double *)sysmem_newptrclear(sizeof(double)*x->NS);
        x->string = (double **)sysmem_newptrclear(sizeof(double *)*x->NS);	// allocating 3 arrays for every string NS
        x->lastString = (double **)sysmem_newptrclear(sizeof(double *)*x->NS);
        x->nextString = (double **)sysmem_newptrclear(sizeof(double *)*x->NS);
            
        x->hammer = (t_pianoHammer **)sysmem_newptrclear(sizeof(t_pianoHammer *)*x->NH);
        x->rattle = (t_rattle **)sysmem_newptrclear(sizeof(t_rattle *)*x->numberOfPreparations);
        x->rubber = (t_rubber **)sysmem_newptrclear(sizeof(t_rubber *)*x->numberOfPreparations);
        x->trap = (t_trap **)sysmem_newptrclear(sizeof(t_trap *)*x->numberOfPreparations);

        for(i=0;i<x->NS;i++)
        {	
            if(x->NS>1)
            {
                x->detuneSpread[i] = x->D*i/(x->NS-1) - x->D/2;
                x->f[i] = x->f0*pow(2, x->detuneSpread[i]/1200.);
            }
            else
                x->f[i] = x->f0;
            
            x->attrDecay = 4;
            x->T30[i] = 4;
            x->attrStiffness = 4;
            x->k[i] = 4;
            x->k_2[i] = pow(x->k[i],2);
            x->attrLoss = 0.002;
            x->b[i] = 0.002;
            x->c[i] = 2*x->f[i];
            x->c_2[i] = pow(x->c[i], 2);
        }
        
        for(i=0;i<x->NS;i++)
        {
            x->string[i] = (double *)sysmem_newptrclear((sizeof(double))*(x->max_N));
            x->lastString[i] = (double *)sysmem_newptrclear((sizeof(double))*(x->max_N));
            x->nextString[i] = (double *)sysmem_newptrclear((sizeof(double))*(x->max_N));
        }
            
        for(i=0;i<x->NH;i++)    //allocating memory for maxnumber of hammers
            x->hammer[i] = t_pianoHammer_new(-0.01, 5000, 1, 10, x->NS);     //initializing hammers.....     dx&dt are set later in dsp-method!!!!
        
        for(i=0;i<x->numberOfPreparations;i++)
        {
            x->rattle[i] = t_rattle_new(3000, 0.9, 0.02, x->NS);             //rattleFundamentalFreq, massDensityRatio, rattleLength
            x->rubber[i] = t_rubber_new(500, 10, 100000, x->NS);
            x->trap[i] = t_trap_new(10000000, 1000, x->NS);
        }
        
        x->rattle1freq = 3000;
        x->rattle1length = 0.02;
        x->rattle1massDensityRatio = 0.9;
    
        x->dt = 1./sys_getsr();        // setting dt and dependent param with "sys_getsr()" provides max from crashing if utility-routines are used before
        x->dt_2 = pow(x->dt,2);		   // turning on processing.....
        
        for(i=0;i<x->NS;i++)
        {
            x->o[i] = (2/x->dt)*(pow(10,(3*x->dt/x->T30[i]))-1);
            x->o_mul_dt_over_2[i] = x->o[i] * x->dt * 0.5;
            x->dxminAll[i] = sqrt(  ((x->c_2[i]*x->dt_2+2*x->b[i]*x->dt) + sqrt( pow(x->c_2[i]*x->dt_2+2*x->b[i]*x->dt,2)+16*x->k_2[i]*x->dt_2)) / 2.);
        }
        
        x->dxmin = t_fmax(x->dxminAll, x->NS);
        x->N = floor(1/x->dxmin);
        x->dx = 1./(double)x->N;
        x->dx_2 = pow(x->dx,2);
        x->dx_4 = pow(x->dx,4);
        
        for(i=0;i<x->NOutlets;i++)
        {
            x->listeningPointNorm[i] = 0.5;
            x->listeningPointIndex[i] = 2+floor(x->listeningPointNorm[i]/x->dx);
        }
        
        for(i=0;i<x->NS;i++)
            calcCoefficients(x, i);
        
        t_pianoHammer_set_dxdt(x->hammer[0], x->dx, x->dt);
        
        for(i=0;i<x->numberOfPreparations;i++)
        {
            
            t_rattle_set_dxdt(x->rattle[i], x->dx, x->dt);
            t_rubber_set_dxdt(x->rubber[i], x->dx, x->dt);
            t_trap_set_dxdt(x->trap[i], x->dx, x->dt);
        }
	}	

	return (x);
}

void ppiano_perform64(ppiano *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	long n = sampleframes;
    long m = 2;
    long i = 0;
     
    double *tmp;
    double *nextString;
    double *string;

    double *stringPtrMinusOne;
    double *stringPtrPlusOne;
    double *stringPtrMinusTwo;
    double *stringPtrPlusTwo;
    double *lastString;
    double *lastStringPtrMinusOne;
    double *lastStringPtrPlusOne;
     
    double *(out[numouts]);
   
    for(i=0;i<numouts;i++)
        out[i] = (double *)(outs[i]);
    
    while (n--) 
    {
        for(i=0;i<x->NS;i++)
        {
            x->string[i][2] = 0; x->string[i][1] = -x->string[i][3]; x->string[i][x->N+2] = 0; x->string[i][x->N+3] = -x->string[i][x->N+1];
            
            nextString = &x->nextString[i][2];
            string = &x->string[i][2];
            lastString =  &x->lastString[i][2];
            stringPtrMinusOne =  &x->string[i][1];
            stringPtrPlusOne = &x->string[i][3];
            stringPtrMinusTwo = &x->string[i][0];
            stringPtrPlusTwo = &x->string[i][4];
            lastStringPtrMinusOne =  &x->lastString[i][1];
            lastStringPtrPlusOne =  &x->lastString[i][3];
            
            for (m=2; m<=x->N; m++)
            {
                *nextString++ = (x->alpha0[i] * *string++) +
                (x->alpha1[i] * (*stringPtrMinusOne++ + *stringPtrPlusOne++)) +
                (x->alpha2[i] * (*stringPtrMinusTwo++ + *stringPtrPlusTwo++)) +
                (x->beta0[i] * *lastString++) +
                (x->beta1[i] * (*lastStringPtrMinusOne++ + *lastStringPtrPlusOne++));
                
            }
        }
        
        
        t_pianoHammer_perform(x->hammer[0], &x->nextString, &x->string);
        
        for(i=0;i<2;i++)
            t_rattle_perform(x->rattle[i], &x->nextString, &x->string);
        
        for(i=0;i<2;i++)
             t_rubber_perform(x->rubber[i], &x->nextString, &x->string);
        
        for(i=0;i<2;i++)
             t_trap_perform(x->trap[i], &x->nextString, &x->string, &x->lastString);
        
        for(i=0;i<numouts;i++)
        {
            *out[i] = 0;
            for(m=0;m<x->NS;m++)
            {
                *out[i] += x->nextString[m][x->listeningPointIndex[i]] * 220; // just some scaling factor which makes sense
                if(x->noteOff)
                    *out[i] = x->lastOut[i] * 0.99999999;
                
                x->lastOut[i] = *out[i];
            }
            
            if(fabs(*out[i]) > 10.)
                resetAll(x);
            
            out[i]++;
        }
 
        for(i=0;i<x->NS;i++)
        {
            tmp = x->lastString[i];
            x->lastString[i] = x->string[i];
            x->string[i] = x->nextString[i];
            x->nextString[i] = tmp;
        }

    }
}

void ppiano_dsp64(ppiano *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    
	int i;
	
	x->dt = 1./samplerate;
	x->dt_2 = pow(x->dt,2);
	
	for(i=0;i<x->NS;i++)
	{
		x->o[i] = (2/x->dt)*(pow(10,(3*x->dt/x->T30[i]))-1);
		x->o_mul_dt_over_2[i] = x->o[i] * x->dt * 0.5;
		x->dxminAll[i] = sqrt(  ((x->c_2[i]*x->dt_2+2*x->b[i]*x->dt) + sqrt( pow(x->c_2[i]*x->dt_2+2*x->b[i]*x->dt,2)+16*x->k_2[i]*x->dt_2)) / 2.);
	}
	
	x->dxmin = t_fmax(x->dxminAll, x->NS);
	x->N = floor(1/x->dxmin);
	x->dx = 1./(double)x->N;
	x->dx_2 = pow(x->dx,2);
	x->dx_4 = pow(x->dx,4);
	
	for(i=0;i<x->NOutlets;i++)
		x->listeningPointIndex[i] = 2+floor(x->listeningPointNorm[i]/x->dx);
    
	for(i=0;i<x->NS;i++)
        calcCoefficients(x, i);
	
    t_pianoHammer_set_dxdt(x->hammer[0], x->dx, x->dt);
    
	for(i=0;i<x->numberOfPreparations;i++)
    {
		t_rattle_set_dxdt(x->rattle[i], x->dx, x->dt);
		t_rubber_set_dxdt(x->rubber[i], x->dx, x->dt);
        t_trap_set_dxdt(x->trap[i], x->dx, x->dt);
    }
	   
    object_method(dsp64, gensym("dsp_add64"), x, ppiano_perform64, 0, NULL);
    
}

//------------------------------------------------------------------------------------------------------------------------------------------//

//																-perform method-

//------------------------------------------------------------------------------------------------------------------------------------------//

void ppiano_free(ppiano *x)
{
	int i;
	for(i=0;i<x->NS;i++)
	{	
		sysmem_freeptr(x->string[i]);
		sysmem_freeptr(x->lastString[i]);
		sysmem_freeptr(x->nextString[i]);
	}
    
	sysmem_freeptr(x->string);
	sysmem_freeptr(x->lastString);
	sysmem_freeptr(x->nextString);
	sysmem_freeptr(x->T30);
	sysmem_freeptr(x->k);				
	sysmem_freeptr(x->k_2);   				
	sysmem_freeptr(x->b);    				
	sysmem_freeptr(x->o);   					
	sysmem_freeptr(x->o_mul_dt_over_2);  	
	sysmem_freeptr(x->f); 
	sysmem_freeptr(x->detuneSpread); 
	sysmem_freeptr(x->c); 
	sysmem_freeptr(x->c_2); 
	sysmem_freeptr(x->dxminAll); 
	sysmem_freeptr(x->alpha0); 
	sysmem_freeptr(x->alpha1); 
	sysmem_freeptr(x->alpha2); 
	sysmem_freeptr(x->beta0); 
	sysmem_freeptr(x->beta1); 
	sysmem_freeptr(x->listeningPointNorm); 
	sysmem_freeptr(x->listeningPointIndex);
	
	for(i=0;i<x->NH;i++)
		t_pianoHammer_free(x->hammer[i]);
	for(i=0;i<x->numberOfPreparations;i++)
    {
		t_rattle_free(x->rattle[i]);
		t_rubber_free(x->rubber[i]);
		t_trap_free(x->trap[i]);
    }
	
	sysmem_freeptr(x->hammer);	
	sysmem_freeptr(x->rattle);
	sysmem_freeptr(x->rubber);
	sysmem_freeptr(x->trap);
}

void ppiano_attr_set_listeningPoint1(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    setListeningPoint(x, 0, atom_getfloat(argv));
}

t_max_err ppiano_attr_get_listeningPoint1(ppiano *x, t_object *attr, long *argc, t_atom ** argv)
{
    char alloc;
    atom_alloc(argc, argv, &alloc);
    if(x->NOutlets >= 1)
        atom_setfloat(*argv, x->listeningPointNorm[0]);
    
    return 0;
}

void ppiano_attr_set_listeningPoint2(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    setListeningPoint(x, 1, atom_getfloat(argv));
}

t_max_err ppiano_attr_get_listeningPoint2(ppiano *x, t_object *attr, long *argc, t_atom ** argv)
{
    char alloc;
    atom_alloc(argc, argv, &alloc);
    if(x->NOutlets >= 2)
        atom_setfloat(*argv, x->listeningPointNorm[1]);
    
    return 0;
}

void ppiano_attr_set_trap1(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->trap1 = atom_getlong(argv);
    trapOn(x, 1, x->trap1);
}

void ppiano_attr_set_trap1position(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->trap1position = atom_getfloat(argv);
    setTrapPosition(x, 1, x->trap1position);
}

void ppiano_attr_set_trap2(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->trap2 = atom_getlong(argv);
    trapOn(x, 2, x->trap2);
}

void ppiano_attr_set_trap2position(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->trap2position  = atom_getfloat(argv);
    setTrapPosition(x, 2, x->trap2position);
}

void ppiano_attr_set_rubber1(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rubber1 = atom_getlong(argv);
    rubberOn(x, 1, x->rubber1 );
}

void ppiano_attr_set_rubber1position(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rubber1position = atom_getfloat(argv);
    setRubberPosition(x, 1, x->rubber1position);
}

void ppiano_attr_set_rubber2(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rubber2  = atom_getlong(argv);
    rubberOn(x, 2, x->rubber2);
}

void ppiano_attr_set_rubber2position(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rubber2position = atom_getfloat(argv);
    setRubberPosition(x, 2, x->rubber2position);
}

void ppiano_attr_set_rattle1(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rattle1 = atom_getlong(argv);
    rattleOn(x, 1, x->rattle1);
}

void ppiano_attr_set_rattle1position(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rattle1position = atom_getfloat(argv);
    setRattlePosition(x, 1, x->rattle1position );
}

void ppiano_attr_set_rattle1freq(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rattle1freq = atom_getfloat(argv);
    setRattleFrequency(x, 1, x->rattle1freq );
}

void ppiano_attr_set_rattle2freq(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rattle2freq = atom_getfloat(argv);
    setRattleFrequency(x, 2, x->rattle2freq );
}

void ppiano_attr_set_rattle1massDensityRatio(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rattle1massDensityRatio = atom_getfloat(argv);
    setRattleMassDensityRatio(x, 1, x->rattle1massDensityRatio);
}

void ppiano_attr_set_rattle1length(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rattle1length= atom_getfloat(argv);
    setRattleLength(x, 1, x->rattle1length);
}

void ppiano_attr_set_rattle2massDensityRatio(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rattle2massDensityRatio = atom_getfloat(argv);
    setRattleMassDensityRatio(x, 2, x->rattle1massDensityRatio);
}

void ppiano_attr_set_rattle2length(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rattle2length= atom_getfloat(argv);
    setRattleLength(x, 2, x->rattle1length);
}

void ppiano_attr_set_rattle2(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rattle2 = atom_getlong(argv);
    rattleOn(x, 2, x->rattle2);
}

void ppiano_attr_set_rattle2position(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->rattle2position  = atom_getfloat(argv);
    setRattlePosition(x, 2, x->rattle2position);
}

void ppiano_attr_set_stringDecay(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->attrDecay  = atom_getfloat(argv);
    setStringDecay(x, x->attrDecay);
}

void ppiano_attr_set_stringStiffness(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->attrStiffness  = atom_getfloat(argv);
    setStringStiffness(x, x->attrStiffness);
}

void ppiano_attr_set_stringLoss(ppiano *x, t_object *attr, short argc, t_atom *argv)
{
    x->attrLoss  = atom_getfloat(argv);
    setStringLoss(x, x->attrLoss);
}

int C74_EXPORT main(void)
{
    post("tr.ppiano~ v0.8");
    post("©2015 ©2016 based on Dr. Stefan Bilbao`s Prepared Piano Sound Synthesis");
    post("Ported from Matlab to Max/MSP by Thomas Resch");
    
    t_class *c;
    
    c = class_new("tr.ppiano~", (method)ppiano_new, (method)ppiano_free, (short)sizeof(ppiano), 0L,
                  A_GIMME, 0);
    class_addmethod(c,(method)ppiano_dsp64,		"dsp64",	A_CANT, 0);		// New 64-bit MSP dsp chain compilation for Max 6
    class_addmethod(c,(method)setHammerPosition, "setHammerPosition", A_DEFLONG, A_DEFFLOAT, 0);
    class_addmethod(c,(method)setNumberOfStrings, "setNumberOfStrings", A_DEFLONG,0);
    class_addmethod(c,(method)setListeningPoint, "setlisteningpoint", A_DEFLONG, A_DEFFLOAT, 0);
    class_addmethod(c,(method)setStringParameters,"setStringParameters", A_GIMME, 0);
    class_addmethod(c,(method)setRubberParameters,"setRubberParameters", A_GIMME, 0);
    class_addmethod(c,(method)setFrequency, "setFrequency", A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(c,(method)setTrapParameters, "setTrapParameters", A_GIMME, 0);
    class_addmethod(c,(method)noteon, "list", A_GIMME, 0);
    class_addmethod(c,(method)noteOff, "noteoff", 0);
    class_addmethod(c,(method)resetAll, "reset", 0);
    class_addmethod(c,(method)plug1, "plug", A_GIMME, 0);
    
    CLASS_ATTR_CHAR(c, "pitchasmidinote", 0, ppiano,  pitchAsMidi);
    CLASS_ATTR_STYLE_LABEL(c, "pitchasmidinote", 0, "onoff", "Pitch as Midi Note");
    CLASS_ATTR_SAVE(c, "pitchasmidinote", 0);
    
    CLASS_ATTR_FLOAT(c, "listeningpoint1", 0, ppiano,  floatDummy);
    CLASS_ATTR_LABEL(c, "listeningpoint1", 0, "Listening Point 1");
    CLASS_ATTR_SAVE(c, "listeningpoint1", 0);
    CLASS_ATTR_ACCESSORS(c, "listeningpoint1", ppiano_attr_get_listeningPoint1, ppiano_attr_set_listeningPoint1);
    
    CLASS_ATTR_FLOAT(c, "listeningpoint2", 0, ppiano,  floatDummy);
    CLASS_ATTR_LABEL(c, "listeningpoint2", 0, "Listening Point 2");
    CLASS_ATTR_SAVE(c, "listeningpoint2", 0);
    CLASS_ATTR_ACCESSORS(c, "listeningpoint2", ppiano_attr_get_listeningPoint2, ppiano_attr_set_listeningPoint2);
    
    CLASS_ATTR_FLOAT(c, "stringdecay", 0, ppiano,  attrDecay);
    CLASS_ATTR_LABEL(c, "stringdecay", 0, "String Decay");
    CLASS_ATTR_SAVE(c, "stringdecay", 0);
    CLASS_ATTR_ACCESSORS(c, "stringdecay", NULL, ppiano_attr_set_stringDecay);
    
    CLASS_ATTR_FLOAT(c, "stringstiffness", 0, ppiano,  attrStiffness);
    CLASS_ATTR_LABEL(c, "stringstiffness", 0, "String Stiffness");
    CLASS_ATTR_SAVE(c, "stringstiffness", 0);
    CLASS_ATTR_ACCESSORS(c, "stringstiffness", NULL, ppiano_attr_set_stringStiffness);
    
    CLASS_ATTR_FLOAT(c, "stringfreqdependentlosses", 0, ppiano,  attrLoss);
    CLASS_ATTR_LABEL(c, "stringfreqdependentlosses", 0, "String Frequency Dependent Losses");
    CLASS_ATTR_SAVE(c, "stringfreqdependentlosses", 0);
    CLASS_ATTR_ACCESSORS(c, "stringfreqdependentlosses", NULL, ppiano_attr_set_stringLoss);
    
    CLASS_ATTR_CHAR(c, "trap1", 0, ppiano,  trap1);
    CLASS_ATTR_STYLE_LABEL(c, "trap1", 0, "onoff", "Trap 1 On/Off");
    CLASS_ATTR_SAVE(c, "trap1", 0);
    CLASS_ATTR_ACCESSORS(c, "trap1", NULL, ppiano_attr_set_trap1);
    
    CLASS_ATTR_FLOAT(c, "trap1position", 0, ppiano,  trap1position);
    CLASS_ATTR_LABEL(c, "trap1position", 0, "Trap 1 Position");
    CLASS_ATTR_SAVE(c, "trap1position", 0);
    CLASS_ATTR_ACCESSORS(c, "trap1position", NULL, ppiano_attr_set_trap1position);
    
    CLASS_ATTR_CHAR(c, "trap2", 0, ppiano,  trap2);
    CLASS_ATTR_STYLE_LABEL(c, "trap2", 0, "onoff", "Trap 2 On/Off");
    CLASS_ATTR_SAVE(c, "trap2", 0);
    CLASS_ATTR_ACCESSORS(c, "trap2", NULL, ppiano_attr_set_trap2);
    
    CLASS_ATTR_FLOAT(c, "trap2position", 0, ppiano,  trap2position);
    CLASS_ATTR_LABEL(c, "trap2position", 0, "Trap 2 Position");
    CLASS_ATTR_SAVE(c, "trap2position", 0);
    CLASS_ATTR_ACCESSORS(c, "trap2position", NULL, ppiano_attr_set_trap2position);
    
    CLASS_ATTR_CHAR(c, "rubber1", 0, ppiano,  rubber1);
    CLASS_ATTR_STYLE_LABEL(c, "rubber1", 0, "onoff", "Rubber 1 On/Off");
    CLASS_ATTR_SAVE(c, "rubber1", 0);
    CLASS_ATTR_ACCESSORS(c, "rubber1", NULL, ppiano_attr_set_rubber1);
    
    CLASS_ATTR_FLOAT(c, "rubber1position", 0, ppiano,  rubber1position);
    CLASS_ATTR_LABEL(c, "rubber1position", 0, "Rubber 1 Position");
    CLASS_ATTR_SAVE(c, "rubber1position", 0);
    CLASS_ATTR_ACCESSORS(c, "rubber1position", NULL, ppiano_attr_set_rubber1position);
    
    CLASS_ATTR_CHAR(c, "rubber2", 0, ppiano,  rubber2);
    CLASS_ATTR_STYLE_LABEL(c, "rubber2", 0, "onoff", "Rubber 2 On/Off");
    CLASS_ATTR_SAVE(c, "rubber2", 0);
    CLASS_ATTR_ACCESSORS(c, "rubber2", NULL, ppiano_attr_set_rubber2);
    
    CLASS_ATTR_FLOAT(c, "rubber2position", 0, ppiano,  rubber2position);
    CLASS_ATTR_LABEL(c, "rubber2position", 0, "Rubber 2 Position");
    CLASS_ATTR_SAVE(c, "rubber2position", 0);
    CLASS_ATTR_ACCESSORS(c, "rubber2position", NULL, ppiano_attr_set_rubber2position);
    
    CLASS_ATTR_CHAR(c, "rattle1", 0, ppiano,  rattle1);
    CLASS_ATTR_STYLE_LABEL(c, "rattle1", 0, "onoff", "Rattle 1 On/Off");
    CLASS_ATTR_SAVE(c, "rattle1", 0);
    CLASS_ATTR_ACCESSORS(c, "rattle1", NULL, ppiano_attr_set_rattle1);
    
    CLASS_ATTR_FLOAT(c, "rattle1position", 0, ppiano,  rattle1position);
    CLASS_ATTR_LABEL(c, "rattle1position", 0, "Rattle 1 Position");
    CLASS_ATTR_SAVE(c, "rattle1position", 0);
    CLASS_ATTR_ACCESSORS(c, "rattle1position", NULL, ppiano_attr_set_rattle1position);
    
    CLASS_ATTR_FLOAT(c, "rattle1fundamentalfreq", 0, ppiano,  rattle1freq);
    CLASS_ATTR_LABEL(c, "rattle1fundamentalfreq", 0, "Rattle 1 Fundamental Frequency");
    CLASS_ATTR_SAVE(c, "rattle1fundamentalfreq", 0);
    CLASS_ATTR_ACCESSORS(c, "rattle1fundamentalfreq", NULL, ppiano_attr_set_rattle1freq);
    
    CLASS_ATTR_FLOAT(c, "rattle1massdensityratio", 0, ppiano,  rattle1massDensityRatio);
    CLASS_ATTR_LABEL(c, "rattle1massdensityratio", 0, "Rattle 1 Mass Density Ratio");
    CLASS_ATTR_SAVE(c, "rattle1massdensityratio", 0);
    CLASS_ATTR_ACCESSORS(c, "rattle1massdensityratio", NULL, ppiano_attr_set_rattle1massDensityRatio);
    
    CLASS_ATTR_FLOAT(c, "rattle1length", 0, ppiano,  rattle1length);
    CLASS_ATTR_LABEL(c, "rattle1length", 0, "Rattle 1 Length");
    CLASS_ATTR_SAVE(c, "rattle1length", 0);
    CLASS_ATTR_ACCESSORS(c, "rattle1length", NULL, ppiano_attr_set_rattle1length);
    
    CLASS_ATTR_FLOAT(c, "rattle2fundamentalfreq", 0, ppiano,  rattle2freq);
    CLASS_ATTR_LABEL(c, "rattle2fundamentalfreq", 0, "Rattle 2 Fundamental Frequency");
    CLASS_ATTR_SAVE(c, "rattle2fundamentalfreq", 0);
    CLASS_ATTR_ACCESSORS(c, "rattle2fundamentalfreq", NULL, ppiano_attr_set_rattle2freq);
    
    CLASS_ATTR_FLOAT(c, "rattle2massdensityratio", 0, ppiano,  rattle2massDensityRatio);
    CLASS_ATTR_LABEL(c, "rattle2massdensityratio", 0, "Rattle 2 Mass Density Ratio");
    CLASS_ATTR_SAVE(c, "rattle2massdensityratio", 0);
    CLASS_ATTR_ACCESSORS(c, "rattle2massdensityratio", NULL, ppiano_attr_set_rattle1massDensityRatio);
    
    CLASS_ATTR_FLOAT(c, "rattle2length", 0, ppiano,  rattle2length);
    CLASS_ATTR_LABEL(c, "rattle2length", 0, "Rattle 2 Length");
    CLASS_ATTR_SAVE(c, "rattle2length", 0);
    CLASS_ATTR_ACCESSORS(c, "rattle2length", NULL, ppiano_attr_set_rattle2length);
    
    CLASS_ATTR_CHAR(c, "rattle2", 0, ppiano,  rattle2);
    CLASS_ATTR_STYLE_LABEL(c, "rattle2", 0, "onoff", "Rattle 2 On/Off");
    CLASS_ATTR_SAVE(c, "rattle2", 0);
    CLASS_ATTR_ACCESSORS(c, "rattle2", NULL, ppiano_attr_set_rattle2);
    
    CLASS_ATTR_FLOAT(c, "rattle2position", 0, ppiano,  rattle2position);
    CLASS_ATTR_LABEL(c, "rattle2position", 0, "Rattle 2 Position");
    CLASS_ATTR_SAVE(c, "rattle2position", 0);
    CLASS_ATTR_ACCESSORS(c, "rattle2position", NULL, ppiano_attr_set_rattle2position);
    
    class_dspinit(c);
    class_register(CLASS_BOX, c);
    ppiano_class = c;
}





