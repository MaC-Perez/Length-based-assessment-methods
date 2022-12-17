//MULTIFAN
//Fournier et al. 1990 CJFAS 47 301-317	
//Schnute & Fournier 1980 CJFAS 37: 1337-1351
//Nesslage 2022, modified from code by Robert Ahrens

DATA_SECTION
	//read in info from dat file
	init_int nsets; //# of lf data sets (Fournier's Na)
	init_int nbins; //# of length bins (Fournier's Nl)
	init_number minbin; //midpt of smallest bin
	init_number binw; //bin width 
	init_ivector years(1,nsets); //years represented by each data set
	init_ivector months(1,nsets); //months represented by each data set
	init_matrix lfdata(1,nsets,1,nbins); //lf data 
	init_int eof; //read-in error check
	int iter;
	!!iter=0;

	//!!cout<<"nsets"<<nsets<<endl;
	//!!cout<<"nbins"<<nbins<<endl;
	//!!cout<<"minbin"<<minbin<<endl;
	//!!cout<<"binw"<<binw<<endl;
	//!!cout<<"years"<<years<<endl;
	//!!cout<<"months"<<months<<endl;
	//!!cout<<"lfdata"<<lfdata<<endl;
	//!!cout<<"eof"<<eof<<endl;
	//!!exit(0);
	
	LOCAL_CALCS
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS

	//switch to reading in info from control file
	!!ad_comm::change_datafile_name("SBFT.ctl"); 
	init_number fage; //first age 
	init_int nage; //number of ages
	init_int vbkph;//phase for vbk
	init_number ivbk; //initial guess for vonB K value
	init_int lfirstph;//phase for lfirst
	init_number ilfirst; // initial guess for mean length for first age
	init_int dlph;//phase for dl
	init_number idl; // initial guess for mean length for last age
	init_int lam1ph;//phase for lam1
	init_number ilam1; // initial guess for lambda 1 (magnitude of sds for meanL@A)
	init_int lam2ph;//phase for lam2
	init_number ilam2; // initial guess for lambda 2 (length-dep trend in sd for meanL@A)
	init_int tauph;//phase for tau	
	init_number itau; // initial guess for tau (var of sampling errors)
	init_int pmeanlonesw; //switch for meanLone
	init_number pmeanlone; //prior mean for first age class
	init_number psdmeanlone;//prior sd for first age class
	init_int pmeanltwosw; //switch for meanLtwo
	init_number pmeanl2;//prior mean for second age class
	init_number psdmeanl2;//prior sd for second age class
	init_int pmeanvbksw; //switch for vbk
	init_number pmeanvbk;//prior mean for vbk
	init_number psdmeanvbk;//prior sd for vbk
	init_int pmeanlinfsw; //switch for meanLtwo
	init_number pmeanlinf;//prior mean for linf
	init_number psdmeanlinf;//prior sd for linf
	init_number multinomform;//choice of multinomial (1=regular, 2=Fournier robust)
	init_int eof2;
	
	
	LOCAL_CALCS
		if(eof2!=999)
		{
			cout<<"Error reading control.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS
	vector age(1,nage);
	vector lbinmids(1,nbins);
	matrix opl(1,nsets,1,nbins);
	LOCAL_CALCS
		age.fill_seqadd(fage,1.);
		lbinmids.fill_seqadd(minbin,binw);
		//opl=lfdata;
		for(int i=1;i<=nsets;i++)opl(i)=lfdata(i)/(rowsum(lfdata)(i));
	END_CALCS

	//!!cout<<"age"<<age<<endl;
	//!!cout<<"lbinmids"<<lbinmids<<endl;
	//!!cout<<"opl"<<opl<<endl;
	//!!exit(0);
	
	
PARAMETER_SECTION
	//init_number vbk(vbkph); // vonB K parameter
	init_bounded_number log_vbk(-5,0,vbkph); // vonB K parameter
	init_number log_lfirst(lfirstph); //mean length for first age
	init_number log_dl(dlph); //mean length for first age
	init_number log_lam1(lam1ph);
	init_number log_lam2(lam2ph);
	init_number log_tau(tauph);
	init_bounded_matrix log_npage(1,nsets,2,nage,-10,10); 
	vector pl(1,nbins);
	!!log_vbk=log(ivbk);
	!!log_lfirst=log(ilfirst);
	!!log_dl=log(idl);
	!!log_lam1=log(ilam1);
	!!log_lam2=log(ilam2);
	!!log_tau=log(itau);
	number vbk;
	number npar;
	number AIC;
	number AICc;
	number BIC
	number lfirst; //mean length first age class
	number lone; //mean length of first age class in month 1
	number ltwo; //mean length at second age class on vonB curve in month 1
	number lmax; // mean length for last age
	//number cv;
	number lam1; //lambda 1 (magnitude of sd)
	number lam2; //lambda 2 (length-dependent trend in the sd; if = 0, sd is indep of length)
	number tau; //parameter defining variance of sampling errors for each data set
	number tzero; //vonB t0 parameter
	number rho; //vonB K parameter =-ln(rho)
	number linf; //vonB Linf parameter
	number fpen; 
	sdreport_number ratio;
	number comp1;
	number comp2;
	number comp3;
	number comp4;
	number comp5;
	number comp6;
	number comp7;
	objective_function_value nll;

	vector la(1,nage);
	vector sdla(1,nage);
	matrix npage(1,nsets,1,nage);
	matrix pla(1,nbins,1,nage);
	matrix ppl(1,nsets,1,nbins);
	matrix pnl(1,nsets,1,nbins)
	matrix meanl(1,nsets,1,nage)
	matrix sdl(1,nsets,1,nage)
	matrix epsilon(1,nsets,1,nbins);

PRELIMINARY_CALCS_SECTION

PROCEDURE_SECTION
	initialization();
	lfset_calculations();
	objective_function();

	if(mceval_phase())
	{
		mcmc_output(); 
	}
	if(last_phase())
	{
	} 

FUNCTION initialization
	fpen=0;
	lmax=mfexp(log_dl)+0.5*binw;//mfexp(log_lfirst)+mfexp(log_dl);
	lfirst=mfexp(log_lfirst);
	lam1=mfexp(log_lam1);
	lam2=mfexp(log_lam2);
	//cv=mfexp(log_cv);
	tau=mfexp(log_tau);
	vbk=mfexp(log_vbk);
	rho=mfexp(-vbk);
	linf=(lmax-lfirst*pow(rho,(nage-fage)))/(1-pow(rho,(nage-fage)));	
	tzero=-(posfun(-(fage-1./log(rho)*log((lmax-lfirst)/(lmax-lfirst*pow(rho,(nage-fage))))),0.001,fpen));
				
FUNCTION lfset_calculations
	dvariable z1;
	dvariable z2;	
	for(int i=1;i<=nsets;i++)
	{
		//npage(i)=page(i)/sum(page(i));//realpage(i);
		npage(i,1)=1.;
		for(int j=2;j<=nage;j++) //loop over ages
		{
		  npage(i,j)=mfexp(log_npage(i,j));
		}
		npage(i)=npage(i)/sum(npage(i));
		la=lfirst+(lmax-lfirst)*(1.-pow(rho,age-fage+(months(i)-1.)/12.))/(1.-pow(rho,nage-fage));
		meanl(i)=la; //Eqn 3.3 - mean length for each age class for each data set
		sdla=lam1*mfexp(lam2*(-1.+2.*(1.-pow(rho,age-fage+(months(i)-1.)/12.))/(1.-pow(rho,nage-fage)))); //Eqn 3.4 - sd of the length distribution of each age class for each data set
		sdl(i)=sdla;
		
		for(int j=1;j<=nage;j++) //loop over ages
		{
			 for(int k=1;k<=nbins;k++) //loop over length bins
			{
				z1=((lbinmids(k)-0.5*binw)-la(j))/sdla(j); 
				z2=((lbinmids(k)+0.5*binw)-la(j))/sdla(j);
				pla(k,j)=cumd_norm(z2)-cumd_norm(z1); 
				//Eqn. 3.2 qijalpha = prob fish of this age lies in interval Z1 to Z2
			}//end nbins
		}//end nage
		
		pl=pla*npage(i); //estimated proportions in each length bin
		ppl(i)=pl/sum(pl);  //standardize so the row sums to one (convert to proportions)
		pnl(i)=ppl(i)*(rowsum(lfdata)(i));  //calculate predicted numbers at length
		epsilon(i)=elem_prod((1-ppl(i)),ppl(i))+0.1/double(nbins);
	}
	lone=linf*(1.-mfexp(-vbk*(age(1)-tzero)));
	ltwo=linf*(1.-mfexp(-vbk*(age(2)-tzero)));
		

FUNCTION objective_function 

	//Multinomial
	if(multinomform==1){	
	comp1=0.;
	comp1=-100.*sum(elem_prod(opl,log(ppl+0.01)));
	comp2=0.;
	comp3=0.;
	comp4=0.;
	if(pmeanlonesw==1)
	{
		comp4=dlnorm(lone,log(pmeanlone),psdmeanlone); //prior on length for first age class;
	}
	comp5=0.;
	if(pmeanltwosw==1)
	{
		comp5=dlnorm(ltwo,log(pmeanl2),psdmeanl2);  //prior on length for first age class
	}
	comp6=0.;
	if(pmeanvbksw==1)
	{
		comp6=dlnorm(vbk,log(pmeanvbk),psdmeanvbk);  //prior for K
	}		
	comp7=0.;
	if(pmeanlinfsw==1)
	{
		comp7=dlnorm(linf,log(pmeanlinf),psdmeanlinf);  //prior for Linf
	}
	}
	
	//Fournier robust multinomial
	if(multinomform==2){	
	
	cout<<"Turn multinomial switch back to 1 and try again!"<<endl;
	exit(0);
	tau=sqrt(1./(nsets*double(nbins))*sum(elem_div(elem_prod(opl-ppl,opl-ppl),elem_prod(ppl,1-ppl))));
	
	comp1=0.;
	comp1=0.5*sum(log(epsilon));
	comp2=0.;
	comp2=double(nbins)*log(tau);
	comp3=0.;
	comp3=sum(log(mfexp(-elem_div(elem_prod((opl-ppl),(opl-ppl)),(2.*epsilon*square(tau)))+0.01)));
	comp4=0.;
	if(pmeanlonesw==1)
	{
		comp4=dlnorm(lone,log(pmeanlone),psdmeanlone); //prior on length for first age class;
	}
	comp5=0.;
	if(pmeanltwosw==1)
	{
		comp5=dlnorm(ltwo,log(pmeanl2),psdmeanl2);  //prior on length for first age class
	}
	comp6=0.;
	if(pmeanvbksw==1)
	{
		comp6=dlnorm(vbk,log(pmeanvbk),psdmeanvbk);  //prior for K
	}		
	comp7=0.;
	if(pmeanlinfsw==1)
	{
		comp7=dlnorm(linf,log(pmeanlinf),psdmeanlinf);  //prior for Linf
	}
	}
	
	nll=comp1+comp2+comp3+comp4+comp5+comp6+comp7+fpen;  //total objective function
	npar=nsets*nage+6.;
	AIC=2*npar-2.*(comp1+comp2+comp3);
	AICc=AIC+(2*npar*(npar-1.))/(nsets*double(nbins)-npar-1.);
	BIC=npar*log(nsets*double(nbins))-2.*(comp1+comp2+comp3);
		//cout<<"comp1"<<comp1<<endl;
		//cout<<"comp2"<<comp2<<endl;
		//cout<<"comp3"<<comp3<<endl;
		//cout<<"comp4"<<comp4<<endl;
		//cout<<"comp5"<<comp5<<endl;
		//cout<<"comp6"<<comp6<<endl;
		//cout<<"comp7"<<comp7<<endl;


//FUNCTION dvariable dlnorm(dvariable param, double log_mean, double sd)
//	return(0.5*square((log(param)-log_mean)/sd));


FUNCTION mcmc_output

	if(iter==0)
	{
		ofstream ofs("par.mcmc");
		ofs<<"vbk\t linf\t to\t"<<endl;
	}
	iter++;
	ofstream ofs("par.mcmc",ios::app);
	ofs<<vbk<<"\t"<<linf<<"\t"<<tzero<<endl;

REPORT_SECTION
  report<<ppl<<endl;
  
  ofstream ofs("other.rep");
	ofs<<opl<<endl;
	ofs<<pnl<<endl;
	ofs<<vbk<<endl;
	ofs<<linf<<endl;
	ofs<<lfirst<<endl;
	ofs<<lmax<<endl;
	

TOP_OF_MAIN_SECTION
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;


