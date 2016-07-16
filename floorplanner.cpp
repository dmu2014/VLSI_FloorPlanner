#include "hfile.h"
#include <ctime>
#include <fstream>
#include <sstream>
#include <iostream>


using namespace fplan;	
using namespace std;
std::string bm_bench;
std::vector<bm::module_t> bm_modules;
std::vector<bm::signal_t> bm_signals;
std::vector<misc::dbl_pair> bm_module_shapes;
std::vector<misc::cstr_t> bm_module_names;
std::vector<misc::dbl_pair> bm_pin_locations;
module_nets_t bm_mnets;
time_t now_g;
double bm_total = 0;

std::string filename(const std::string &pathname)
{
	std::string ret;
	for (size_t i = pathname.length(); i != 0; --i)
		if ((pathname[i-1] == '\\') || (pathname[i-1] == '/'))
			break;
		else
			ret = std::string(1, pathname[i-1])+ret;
	return ret;
}

bool load_bm()
{
	int n = bm::read_bm(bm_bench, bm_modules, bm_signals);

	bm_bench = filename(bm_bench);

	bm_info(n,bm_modules, bm_signals,bm_module_shapes, bm_module_names,bm_pin_locations,bm_mnets);

	bm_total = 0;
	for (size_t i = 0; i < bm_module_shapes.size(); ++i)
		bm_total += bm_module_shapes[i].first*bm_module_shapes[i].second;

	return true;
}



int main(int argc, char *argv[])
{
	
	time_t now = time(0);
	now_g=now;
	string term; 
	string termnet; 
	string modnet;
	string modidnet; 
	string modid; 
	string width; 
	string height;
	string ul; 
	string ur;
	string netline;
	ifstream inFileterm;
	ifstream inFilenets;
	ifstream inFilemodinfo;
	int found;
	int found2;
	int found3;
	int found4;

	string modline;
	double xl = 0;
	double yl = 0;
	
	double x1 ;
	double y1;
	double x2 ;
	double y2;
	double x3 ;
	double y3 ;
	double x4 ;
	double y4 ;
	double x12 ;
	double y12;
	double x22 ;
	double x32 ;
	double x42 ;
	double y42 ;
	double x52 ;
	double y52;
	double x62 ;
	double x72;
	double x82 ;
	double y82 ;
	double x1e;
	double x2e;
	double x3e;
	double x4e;
	double x5e;
	double x6e;
	double centerx;
	double centery;
	bm_bench = argv[1];
	string T_stop;
	string bname;
	
	if (argv[2] ==NULL)
	{
		cout << "Error in Usage.\nType ./floorplanner BenchmarkName Stop_Temperature\nExample Usage: ./floorplanner B10 0.1 \n " << endl;
		cout << "Press any key to Exit"<<endl;

		getchar();
		return 0;
	}
	T_stop = argv[2];
	bname=argv[1];

	//Read Benchmark Files
	int nb = bm::read_bm(bm_bench, bm_modules, bm_signals);

	bm_bench = filename(bm_bench);

	bm_info(nb,bm_modules, bm_signals,bm_module_shapes, bm_module_names,bm_pin_locations,bm_mnets);

	bm_total = 0;
	for (size_t i = 0; i < bm_module_shapes.size(); ++i)
		bm_total += bm_module_shapes[i].first*bm_module_shapes[i].second;

	//Initialise Pseudo Random Number Generator for reproducible results
	srand(0); 

	//Get ACG Floor plan
	acg_fp_t acg_fp(bm_module_shapes.size(), &bm_module_shapes[0]);
	
	const size_t n = bm_module_shapes.size();
	const size_t m = bm_pin_locations.size();
	std::string oname;

	//Define cost function
	cost_ptr wlHP = acg_fp.get_fplan_interface()->HPWL_cost(bm_mnets, m, &bm_pin_locations[0]);
	cost_ptr area = acg_fp.get_fplan_interface()->area_cost();
	cost_ptr cost = 0.8*area + 0.2*wlHP ;

	if (bname.compare("B10") == 0)
		cost = 0.2*area + 0.8*wlHP ;
	if (bname.compare("B30") == 0)
		cost = 0.8*area + 0.2*wlHP ;
	if (bname.compare("B50") == 0)
		cost = 0.4*area + 0.6*wlHP ;
	if (bname.compare("B100") == 0)
		cost = 0.1*area + 0.9*wlHP ;
	if (bname.compare("B200") == 0)
		cost = 0.9*area + 0.1*wlHP ;
	if (bname.compare("B300") == 0)
		cost = 0.6*area + 0.4*wlHP ;
	
	//Print initial Randomised arrrangement
	fprintf(stdout, "----Initial----\narea %e, deadspace %e\n Annealing in progress...\n",
		area->cost(), (area->cost()-bm_total)/area->cost()*100);

	//Perform Simulated Annealing
	anneal(cost, acg_fp.create_perturbation(0, 0.5), int(n), stdout, 10, 0.90, atof(T_stop.c_str()), 0.90);
	
	misc::dbl_pair shape = acg_fp.get_fplan_interface()->get_shape(); // eval
	double w = shape.first;
	double h = shape.second;
	fprintf(stdout,
		"----Final----\narea %e, deadspace %e, width %f, height%f\n",
		area->cost(), (area->cost()-bm_total)/area->cost()*100, w, h);

	
	//Write Results of Modules
	acg_fp.dump(bm_bench+".pl", &bm_module_names[0]);
	
	
	//Compute Terminal Positions
	FILE *fpterm = fopen("terminals.txt", "w");
	//Create Terminal File
	for (size_t i = size_t(nb); i < bm_modules.size(); ++i)
		fprintf(fpterm, "%s\n",bm_modules[i].name.c_str()); 
	fprintf(fpterm, "\n");
	fclose(fpterm);

	//Create Nets File
	FILE *fpnets = fopen("nets.txt", "w");
	for (size_t i = 0; i < bm_signals.size(); ++i)
	{
		const bm::signal_t &s = bm_signals[i];
		
		
		for (size_t j = 0; j < s.pins.size(); ++j)
			fprintf(fpnets, "%s\t", s.pins[j]->name.c_str());
		fprintf(fpnets, "\n");
	}
	fclose(fpnets);
	double xh = w;
	double yh = h;
	x1=2;
	y1=(h/2) - 2;
	x2=0;
	y2=(h/2) +2;
	x3=(w/2)+(w/8) +2;
	y3=(h/2) - 2;
	x4=(w/2)+2;
	y4=(h/2)+2;
	x12=2;
	y12=2;
	x22=(w/4)+2;
	x32=(w/2) +2;
	x42=(3*w/4)+2;
	y42=2;
	x52 = (3*w/4) +2;
	y52 = (h/2) + 2;
	x62 = (w/2) +2;
	x72 = (w/4) + 2;
	x82 = 2;
	y82 = (h/2) + 2;
	x1e = (w/2)+(w/8) - 2;
	x2e = 2;
	x3e = 2;
	x4e = (w/2) - 2;
	x5e = w-2;
	x6e = w-2;
	double y1ee = 2;
	double y4ee = h -2;
	double y5ee = 2;
	double y7ee = h-2;



	bool flag = false;
	inFileterm.open("terminals.txt");
	inFilenets.open("nets.txt");
	inFilemodinfo.open("moduleinfo.txt");
	std::string b = ".pl";
	std::string bench = argv[1]+b;
	FILE *fpres = fopen(bench.c_str(), "a+");

iter:
while ( inFileterm >> term )
{
  while(getline (inFilenets,netline,'\n'))
  {		
	  
	  found=netline.find(term);
	  if(found==std::string::npos)
		  continue;
	  
	  if (found!=std::string::npos)
	  {
		  found2 = netline.find("sb");
		  if (found2!=std::string::npos)
		  {	
			  string nnl = netline.substr(found2+2);
			  
			  found3 = nnl.find("\t");
			  
			  modidnet = nnl.substr(0,found3);
			  while(inFilemodinfo >> modid >> width >> height >> ul >> ur)
			  {
				  found4 = modid.find(modidnet);
				  if (found4==std::string::npos)
					  continue;
				  if (found4!=std::string::npos)
				  {
					  centerx = atof(ul.c_str()) - atof(width.c_str())/2 ;  
					 
					  
					  centery = atof(ur.c_str()) - atof(height.c_str())/2 ;
					  
					  
					  if(centerx <= w/2 & centery <=h/2)
					  {	
						  if(h*centerx <= w*centery)
						  {
							  if(x1e<=0)
							  {
								  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xl, y1ee);
								  y1ee+=2;
							  }
							  else
							  {
								  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x1e, yl);
								  x1e-=2;
							  }
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
								  
								  
							  
						  }

						  else
						  {
							  
								  if(y1<=0)
									  {
										  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x2e, yl);
										  x2e+=2;
									  }
								  else{
								  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xl, y1);
									  y1-=2;
								  }
							  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
								  
								  
						  }
						  
						
					  }

					  if(centerx <= w/2 & centery >h/2)
					  {		
						  if(centery*w+h*centerx <= h*w)
						  {
									 
										  if (y2>=h)
									  {
										  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x3e, yh);
										  x3e+=2;
									  }
										  else
										  {
											  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xl, y2);
											  y2+=2;
										  }
									 
									  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
								  
						  }

						  else
						  {
										  if(x4e<=0)
										  {
											  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xl, y4ee);
											  y4ee-=2;
										  }
										  else
										  {
										  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x4e, yh);
									      x4e-=2;
										  }
									  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
								  
								  
						  }

					  }
					  if(centerx > w/2 & centery <=h/2)
					  {		
						  if(centery*w+h*centerx <= h*w)
						  {
									  if(x3>=w)
									  {
										  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xh, y5ee);
										  y5ee+=2;
									  }
									  else
									  {
								      fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x3, yl);
									  x3+=2;
									  }
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
								  
						  }

						  else
						  {
									  if(y3<=0)
									  {
										   fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x5e, yl);
										   x5e-=2;
										  
									  }
									  else
									  {
									  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xh, y3);
									  y3-=2;	
									  }
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
								  
								  
						  }
						  

					  }
					   if(centerx > w/2 & centery > h/2)
					  {		
						  if(h*centerx <= w*centery)
						  {
									  if(x4>=w)
									  {
										  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xh, y7ee);
										  y7ee-=2;
									  }
									  else
									  {
								      fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x4, yh);
									  x4+=2;
									  }
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
								  
						  }

						  else
						  {
									  if(y4>=h)
									  {
										  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x6e, yh);
										  x6e-=2;
										  
									  }
									  else
									  {
									  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xh, y4);
									  y4+=2;
									  }
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
								  
						  }
						  

					  }
					   //}
					/*else
						{*/
							/*if(centerx <= w/4 & centery <=h/2)
							{
								if (x12 >= w/4)
									  {
										  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xl, y12);
										  y12+=2;
									  }
								else{
								fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x12, yl);
									  x12+=2;	
								}
									  									  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
							}

							if(centerx <= w/4 & centery >h/2)
							{
								if (x82 >= w/4)
									  {
										  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xl, y82);
										  y82+=2;
									  }
								else{
								fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x82, yh);
									  x82+=2;	
								}
									  									  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
							}

							if(centerx <= w/2 & centery <=h/2)
							{
								
								fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x22, yl);
									  x22+=2;	
									  									  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
							}
							if(centerx <= w/2 & centery > h/2)
							{
								
								fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x72, yh);
									  x72+=2;	
									  									  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
							}
							if(centerx <= 3*w/4 & centery <=h/2)
							{
								
								fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x32, yl);
									  x32+=2;	
									  									  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
							}
							if(centerx <= 3*w/4 & centery > h/2)
							{
								
								fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x62, yh);
									  x62+=2;	
									  									  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
							}

							if(centerx <= w & centery <=h/2)
							{
								if (x42 >= w)
									  {
										  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xh, y42);
										  y42+=2;
									  }
								else{
								fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x42, yl);
									  x42+=2;	
								}
									  									  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
							}
							if(centerx <= w & centery > h/2)
							{
								if (x52 >= w)
									  {
										  fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), xh, y52);
										  y52+=2;
									  }
								else{
								fprintf(fpres,"%s\t%f\t%f\n",term.c_str(), x52, yh);
									  x52+=2;	
								}
									  									  
									  inFilemodinfo.seekg (0, ios::beg);
									  inFilenets.seekg (0, ios::beg);
									  goto iter;
							}*/

							



						/*}*/







				  }

			  }
			  
			 
		  }
	  }
  }
  
}

	inFileterm.close();
	inFilenets.close();
	inFilemodinfo.close();
	fclose(fpres);
	
	
	time_t then = time(0);
	time_t diff = then - now;
	cout << "Running time:" << diff << "seconds\n" << endl;

	remove( "terminals.txt" );
	remove( "nets.txt" );
	remove( "moduleinfo.txt" );
	cout << "Press any key to Exit"<<endl;

	getchar();
	return 0;
}

