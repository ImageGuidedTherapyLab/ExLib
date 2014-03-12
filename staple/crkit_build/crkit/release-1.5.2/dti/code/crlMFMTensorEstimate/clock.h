
#ifndef h_clck
#define h_clck

#include <time.h>

class Clock
{
	clock_t start_t, end_t;
	time_t real_start_t, real_end_t ;

public:
	void reset() { end_t=0; start_t=0; real_start_t=0; real_end_t=0;}
	void start() 
	{ 
		end_t=0; start_t = clock(); 
		real_end_t = real_start_t = time(NULL); 
	}
	void stop() 
	{ 
		end_t = clock(); 
		real_end_t = time(NULL);
	}


	float laps() 
	{ 

		clock_t cur_end_t= clock();
		float delay= (float)(cur_end_t - start_t);
		delay/=CLOCKS_PER_SEC;//unix               =(delay/=CLK_TCK);     // Borland , windows
		return(delay); 

	}  //en seconde

	float real_laps() 
	{ 
		float delay= (float)difftime(real_end_t, real_start_t);
		return(delay); 
	}  //en seconde


	void formatTime(float seconds, char *szBuffer)
	{
		int h, mn, s;

		h = (int)(seconds/3600.0f);  seconds -= 3600.0f*h;
		mn = (int)(seconds/60.0f);   seconds -= 60.0f*mn;
		s = (int)(seconds+0.5f); 

		if ( h>0 )
			sprintf(szBuffer,"%02dh%02dm%02ds", h, mn, s);
		else if (mn>0 )
			sprintf(szBuffer,"%02dm%02ds", mn, s);
		else
			sprintf(szBuffer,"%02ds", s);
	}


	void cout_laps() 
	{ 
		float delay=0.0, h_f=0.0,mn_f=0.0;
		int h=0,mn=0,s=0;

		delay=laps();
		h_f=delay/3600.0f; h=(int)h_f;    delay-= 3600.0f*h;
		mn_f=delay/60.0f;  mn=(int)mn_f;  delay-= 60.0f*mn;
		s=(int)(delay+0.5); 

		if (h>0) cout<<h<<" h "<<mn<<" mn "<<s<<" s "<<endl;
		else if (mn>0) cout<<mn<<" mn "<<s<<" s "<<endl;
		else cout<<s<<" s "<<endl;
	}


	void print_laps(ofstream &f) 
	{ 
		float delay=0.0, h_f=0.0,mn_f=0.0;
		int h=0,mn=0,s=0;

		delay=laps();
		h_f=delay/3600.0f; h=(int)h_f;    delay-= 3600.0f*h;
		mn_f=delay/60.0f;  mn=(int)mn_f;  delay-= 60.0f*mn;
		s=(int)(delay+0.5f); 

		if (h>0) f<<h<<" h "<<mn<<" mn "<<s<<" s "<<endl;
		else if (mn>0) f<<mn<<" mn "<<s<<" s "<<endl;
		else f<<s<<" s "<<endl;

	}

	Clock():start_t(0), end_t(0){}
};

#endif


