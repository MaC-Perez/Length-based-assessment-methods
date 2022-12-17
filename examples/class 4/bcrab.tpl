//North Carolina catch survey model
//M. Wilberg
//Colton et al. Fish. Oceanogr. 23:2, 132–146

TOP_OF_MAIN_SECTION
 //increase number of estimated parameters
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(200040);
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
 arrmblsize = 10000000;

DATA_SECTION
 //READ IN DATA HERE
 init_int fyear //first year of model
 init_int lyear //last year of model
 //Catch Data
 init_int ftcyear //first year of total catch
 init_int ltcyear // last year of total catch
 init_vector com_TC_obs(ftcyear,ltcyear) //total catch
 //North Carolina Survey Data
 init_int fncyear //first year of NC survey
 init_int lncyear //last year of NC survey
 init_vector Ia_NCobs(fncyear,lncyear) //adult NC
 init_vector Ir_NCobs(fncyear,lncyear) //juv. NC
 init_number Ia_NCsd
 init_number Ir_NCsd
 //Natural Mortality
 init_number set_M 
 //proportion of recreational harvest
 init_number p_rec
 init_int test //EOF number
 //Total Harvest including rec
 init_vector TC_obs(ftcyear,ltcyear) //total catch
 int y //looping variable for year
 

 LOCAL_CALCS
  if (test!=12345)
  {
  cout<<"Data not reading properly"<<endl;
  cout<<"fyear,lyear:"<<fyear<<","<<lyear<<endl;
  cout<<"ftcyear,ltcyear:"<<ftcyear<<","<<ltcyear<<endl;
  cout<<"Total Catch"<<endl<<TC_obs<<endl;
  cout<<"fncyear,lncyear:"<<fncyear<<","<<lncyear<<endl;
  cout<<"Adult NC survey indices"<<endl<<Ia_NCobs<<endl;
  cout<<"Juv NC survey indices"<<endl<<Ir_NCobs<<endl;
  cout<<"M:"<<set_M<<endl;
  cout<<"EOF test:"<<test<<endl;
  exit(1);
  }
  //Calculate Total Catch
  TC_obs=com_TC_obs*(1.+p_rec); //commercial + rec catch
 END_CALCS

PARAMETER_SECTION
init_bounded_number log_N(0.,20.,1) //log initial adult abundance
init_bounded_number log_mean_R(0.,20.,1) //log mean recruitment
init_bounded_dev_vector log_R_devs(fyear,lyear,-20.,20.,1) //log recruitment deviations
number pen //penalty function for N to avoid zero/negative abundance
number qa_NC //Age 1+ catchabilty
number qr_NC //age-0 catchability
vector N(fyear,lyear)
vector R(fyear,lyear)
vector N_NC(fncyear,lncyear)
vector R_NC(fncyear,lncyear)
vector Ia_NCest(fncyear,lncyear) //index of NC adults
vector Ir_NCest(fncyear,lncyear) //index of NC recruits
vector u(fyear,lyear)
number ubar //mean u
//Likelihoods
number La_NC //likelihood for adult index
number Lr_NC //likelihood for recruit index
number u_prior //prior on u
//number Lc
objective_function_value NegLL //total likelihood (sum of all likelihood components)

 LOCAL_CALCS
  //Set initial parameter values
  log_N=5.;
  log_mean_R=6.;
 END_CALCS

PROCEDURE_SECTION
 //fill in recruitment
 R=exp(log_mean_R+log_R_devs);
 //Fill in first year of adult abundance
 N(fyear)=exp(log_N);
 //Fill in rest of years abundance
 pen=.0;
 for(y=fyear;y<lyear;y++)
 {
 N(y+1)=((N(y)+R(y))*exp(-(set_M/2.))-TC_obs(y))*exp(-(set_M/2.));
 if (N(y+1)<=.0)
 {
 pen+=square(N(y+1));
 N(y+1)=1;
 }
 }
 //calculate u
 for(y=fyear;y<=lyear;y++)
 {
 u(y)=TC_obs(y)/((N(y)+R(y))*exp(-set_M/2.));
 }
 ubar=sum(u)/double(lyear-fyear+1);
 //incorporate timing of the surveys
 N_NC=(N);
 R_NC=(R);
 //calculate catchability
 qa_NC=exp(sum(log(Ia_NCobs)-log(N_NC))/double(lncyear-fncyear+1));
 qr_NC=exp(sum(log(Ir_NCobs)-log(R_NC))/double(lncyear-fncyear+1));
 //calculate indices
 Ia_NCest=qa_NC*N_NC;
 Ir_NCest=qr_NC*R_NC;
 //lognormal likelihood for indices of abundance
 La_NC=double(lncyear-fncyear+1)*log(Ia_NCsd)+0.5*norm2(log(Ia_NCobs)-log(Ia_NCest))/square(Ia_NCsd);
 Lr_NC=double(lncyear-fncyear+1)*log(Ir_NCsd)+0.5*norm2(log(Ir_NCobs)-log(Ir_NCest))/square(Ir_NCsd);
 //implement beta distribution prior on u
 //prior=-(alpha-1)*log(u)-(beta-1)*log(1-u)
 u_prior=-(4.-1.)*log(ubar)-(4.-1.)*log(1.-ubar);
 //calculate overall negative log likelihood
 NegLL=La_NC+Lr_NC+pen+u_prior;

REPORT_SECTION
 /*report << "Beginning of report section" << endl;
 report << "Likelihood components" << endl;
 report << "La_NC" << endl;
 report << La_NC << endl;
 report << "Lr_NC" << endl << Lr_NC <<endl;
 report << "penalty" << " " << pen << endl;
 report << "u_prior" << " " << u_prior << endl;
 report << "NEGLL" << " " << NegLL << endl;
 report << endl;
 report << "CVs for indices of abundance" << endl;
 report <<Ia_NCsd << " " << Ir_NCsd << endl;
 report << "M" << endl;
 report << set_M << endl;
 report << "Adult catchability" << endl;
 report << qa_NC<< endl;
 report << "recruit catchability" << endl;
 report << qr_NC << endl;
 report << "Year Adult_N Rec_N u Catch" << endl;
 */
 
 for(y=fyear;y<=lyear;y++)
 {
 report << y << " " << N(y) << " " << R(y) << " " << u(y) << " " << log_R_devs(y) << endl;
 }
 //report << "Year NCadult_est NCadult_obs NCjuv_est NCjuv_obs" << endl;
 
 for(y=fncyear;y<=lncyear;y++)
 {
 report << y << " " << Ia_NCest(y) << " " << Ia_NCobs(y) << " " << Ir_NCest(y) << " " << Ir_NCobs(y) << endl;
 }

RUNTIME_SECTION
 maximum_function_evaluations 5000, 25000, 20000, 20000, 20000, 20000
 //leave empty line below here
 