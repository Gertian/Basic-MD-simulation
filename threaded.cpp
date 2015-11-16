#include <vector>
#include <iostream>
#include <cmath>
#include <math.h> // M_PI
#include <stdlib.h>
#include <fstream>
#include <string>
#include <omp.h> //threading
#include <time.h> //timing
#include <chrono>

using namespace std;

//******************************************************
//
//overloading on the vector class
//
//******************************************************

//overloaded <<
//=============
template<typename T> 
ostream& operator <<(ostream& os, const vector<T>& a){
	unsigned int size = a.size();
	for(unsigned int i = 0; i < size; ++i){
		os << a[i] << "	";
	}
	return os;
}

//overloaded *
//============
template<typename T>
vector<T> operator *(const double a, vector<T> b){
	unsigned int size = b.size();
	for(unsigned int i =0; i<size; ++i){
		b[i] = b[i]*a;
	}
	return b;
}

//overloaded +
//============
template<typename T>
vector<T> operator +(const vector<T> a, vector<T> b){
	unsigned int size = b.size();
	for(unsigned int i = 0; i< size; ++i){
		b[i] = b[i]+a[i];
	}
	return b;
}

//**************************************************************************
//**************************************************************************
//
//
//class vec that represents a single particle with and x, y and z component
//
//
//**************************************************************************
//**************************************************************************
class Vec{
	private:
		double _x,_y,_z;
	public:		
		Vec(){_x=0.; _y=0.; _z=0.;}
		Vec(double x, double y, double z){_x=x;_y=y;_z=z;}
		Vec(const Vec& a){_x=a._x;_y=a._y;_z=a._z;}
	public:
		double getX()const{return _x;}
		double getY()const{return _y;}
		double getZ()const{return _z;}
		double norm()const{return sqrt(_x*_x+_y*_y+_z*_z);}
};

//*************************
//
//OPERERATOR OVERLOADING
//
//*************************

//overloaded +
//============
Vec operator +(Vec a,const Vec b){
	a = Vec(a.getX()+b.getX(),a.getY()+b.getY(),a.getZ()+b.getZ());
	return a;
}

//overloaded *
//============
Vec operator *(Vec a, const double b){
	return Vec(a.getX()*b,a.getY()*b,a.getZ()*b);
}

//overloaded *
//============
double operator *(const Vec a, const Vec b){
	return a.getX()*b.getX()+a.getY()*b.getY()+a.getZ()*b.getZ();
}

//overloaded <<
//=============
ostream& operator <<(ostream& os,const Vec& a){
	os << a.getX() << "	" << a.getY() << "	" << a.getZ();
	return os;
}

//*************************************************************************
//*************************************************************************
//
//
//class Object that holds the position, speed and pair list of objects
//
//
//*************************************************************************
//*************************************************************************
class Object{
	private:
		Vec _position,_speed;
		vector<Object*> _friends;
	public:
		Object(){_position = Vec(); _speed = Vec();}
		Object(Vec position, Vec speed){_position = position; _speed = speed;}
		Object(const Object& a){_position = a._position; _speed = a._speed; _friends = a._friends;}
	public:
		Vec getPosition()const{return _position;}
		Vec getSpeed()const{return _speed;}
		unsigned int getAmountFriends()const{return _friends.size();}
		Object getFriendAt(unsigned int i)const{return *(_friends[i]);}
		vector<Object> getFriends(){
			vector<Object> output;
			for(unsigned int i = 0; i<_friends.size(); ++i){
				output.push_back(*(_friends[i]));
			}
			return output;
		} 
		
		double speedNorm()const{return _speed.norm();}
		double positionNorm()const{return _position.norm();}

		void shiftPosition(Vec a,double systemsize){
			//first shift the position
			_position = _position+a;
			//if the new positions are out of system, shift them back !
			double corrx;
			double corry;
			double corrz;
			//we hold a bool that is true after the while loop only when non of the positions was out of bounds during this iteration
			bool isOK = false;
			while(!isOK){
				//true to put the bool on OK.
				isOK = true;
				corrx = 0;
				corry = 0;
				corrz = 0;
				if(_position.getX() > systemsize){
					corrx = -systemsize;
					//we had to make an adaptation. Make the bool false!
					isOK = false;
				}
				if(_position.getX() < 0){
					corrx = systemsize;
					isOK = false;
				}
				if(_position.getY() > systemsize){
					corry = -systemsize;
					isOK = false;
				}
				if(_position.getY() < 0){
					corry = systemsize;
					isOK = false;
				}
				if(_position.getZ() > systemsize){
					corrz = -systemsize;
					isOK = false;
				}
				if(_position.getZ() < 0){
					corrz = systemsize;
					isOK = false;
				}
				//make the correction in case the object was shifted out of bounds
				_position = _position + Vec(corrx,corry,corrz);
			}

		}
		void shiftSpeed(Vec a){_speed = _speed+a;}
		void multiplySpeed(double a){_speed = _speed*a;}

		void clearFriends(){_friends.clear();}
		void addFriend(Object* a){_friends.push_back(a);}
};

//***********************
//
//OPERERATOR OVERLOADING//
//
//**********************

//overloaded <<
//============
ostream& operator <<(ostream& os,const Object& a){
	os << a.getPosition() << "	" << a.getSpeed();
	return os;
}

//********************************************
//A function that generates two exponential distributed numbers
//********************************************
vector<double> expRand(){
	double rand1 = ((double) rand() / (RAND_MAX));
	double rand2 = ((double) rand() / (RAND_MAX));

	double a = sqrt(-2.*log(rand1)); //-->gaussian distribution
	double phi = 2.*M_PI*rand2;      //-->uniform angular distribution

	double v1 = a*cos(phi); 
	double v2 = a*sin(phi);
	

	vector<double> output;
	output.push_back(v1);
	output.push_back(v2);
		
	return (output);
}

//*****************************************
//Given two coordinates, if one want to determine the distance in between them in the minimum image convention one has to do the following
//1)calculate the actual distance, no fancy tricks
//2)if dist > systemsize/2 than we need to subtract systemsize from this value
//3)if dist < -systemsize/2 than we need to add systemsize to this value
//****************************************
double calcDistance(double coord1, double coord2, double systemsize){
	double distance = coord1 - coord2;
	while(distance > systemsize/2.){
		distance = distance - systemsize;
	}
	while(distance < -systemsize/2.){
		distance = distance + systemsize;
	}
	return distance;
}

//*************************************************************************
//*************************************************************************
//
//
//Class Simulation that holds a vector with all the objects, the time, the amount of particles, and the amount of unit cel
//
//
//*************************************************************************
//*************************************************************************
class Simulation{
	private:
		double _time, _systemsize, _celsize;
		double _density;
		double _rcutof;
		double _rcutofextra;
		unsigned int _M;
		unsigned int _amountobjects;
		unsigned int _subs;
		double _subslength;
		vector<Object> _objects;
		vector<Vec> _initialspeeds;
	public:
		double getTime()const{return _time;}
		double getSystemsize()const{return _systemsize;}
		double getCelsize()const{return _celsize;}
		double getDensity()const{return _density;}
		double getRcutof()const{return _rcutof;}
		double getRcutofextra()const{return _rcutofextra;}
		unsigned int getM()const{return _M;}
		unsigned int getAmountobjects()const{return _amountobjects;}
		unsigned int getSubs()const{return _subs;}
		double getSubslength()const{return _subslength;}
		vector<Object> getObjects()const{return(_objects);}
		Object getObjectAt(unsigned int i)const{return(_objects[i]);}
	public:
		//************************
		//function that returns the vector<Vec> of speeds and positions
		//************************
		vector<Vec> getPositions(){
			vector<Vec> output;
			for(unsigned int i = 0; i<_amountobjects; ++i){
				output.push_back(_objects[i].getPosition());
			}
			return output;
		}

		vector<Vec> getSpeeds(){
			vector<Vec> output;
			for(unsigned int i = 0; i<_amountobjects; ++i){
				output.push_back(_objects[i].getSpeed());
			}
			return output;
		}
	
	public:
		//************************
		//fucntion that empties all the friend lists in the system
		//***********************
		void clearSimulationFriends(){
			for(unsigned int i =0; i< _amountobjects; ++i){
				_objects[i].clearFriends();
			}
		}	
		
	public:
		//***************************
		//initialse the _initialspeeds
		//**************************
		void set_Initialspeeds(){
			_initialspeeds = this->getSpeeds();
		}

	public:	
		//************************
		//functions that shift the positions and speeds of all the objects in the simulation
		//************************
		void shiftPositions(vector<Vec> a){
			for(unsigned int i = 0; i<_amountobjects; ++i){
				_objects[i].shiftPosition(a[i],_systemsize);
			}
		}
		void shiftSpeeds(vector<Vec> a){
			for(unsigned int i = 0; i<_amountobjects; ++i){
				_objects[i].shiftSpeed(a[i]);
			}
		}
		void multiplySpeeds(double a){
			for(unsigned int i = 0; i<_amountobjects; ++i){
				_objects[i].multiplySpeed(a);
			}
		}
		void shiftTime(double t){
			_time += t;
		}


	public:
		//******************************************************************************
		//the initialisation procedure of the simultion is as follows
		//initialze the trivial things like time, systemsize,...	
		//->loop over all the cells in the ocject i[0tocels],j[0tocels],k...
		//	->call the expRand fucntion 3 times and store the 6 initial velocities
		//	->initialize the first two particles  using the first 6 velocities
		// 	->reapeat for the two resting particles 
		//*****************************************************************************
		Simulation(double density, unsigned int M, double rcutof, double rcutofextra, unsigned int subs){
			_time = 0;
			_systemsize = pow(4*M*M*M/density,1./3.);
			_density = density;
			_celsize = _systemsize/M;
			_rcutof = rcutof;
			_rcutofextra = rcutofextra;
			_M = M;
			_subs = subs;
			_subslength = _systemsize/subs;
			_amountobjects = 4*M*M*M;
			for(unsigned int i=0; i<_M; ++i){
				for(unsigned int j=0; j<_M; ++j){
					for(unsigned int k=0; k<_M; k++){
						vector<double> speeds;
						//call the exprand function 6 times in a loop;
						for(unsigned int counter=0; counter < 6; counter++){
							vector<double> temp = expRand();
							speeds.push_back(temp.at(0));
							speeds.push_back(temp.at(1));
						}
						//now we have calculated all the needed speeds for the initialisation of this cel.
						//now create the 4 particles in the cel keeping a little ofset into account!
						Object a = Object(Vec(_celsize*(i+1/8.    ),_celsize*(j+1/8.    ),_celsize*(k+1/8.    )),Vec(speeds.at(0),speeds.at(1 ),speeds.at(2 )));
	
						Object b = Object(Vec(_celsize*(i+1/8.+0.5),_celsize*(j+1/8.+0.5),_celsize*(k+1/8.    )),Vec(speeds.at(3),speeds.at(4 ),speeds.at(5 )));					

						Object c = Object(Vec(_celsize*(i+1/8.+0.5),_celsize*(j+1/8.    ),_celsize*(k+1/8.+0.5)),Vec(speeds.at(6),speeds.at(7 ),speeds.at(8 )));					

						Object d = Object(Vec(_celsize*(i+1/8.    ),_celsize*(j+1/8.+0.5),_celsize*(k+1/8.+0.5)),Vec(speeds.at(9),speeds.at(10),speeds.at(11)));					
						_objects.push_back(a);
						_objects.push_back(b);
						_objects.push_back(c);
						_objects.push_back(d);					
					}
				}
			}
			_initialspeeds =  vector<Vec>(_amountobjects,Vec());	
		}
		
	public:
		//******************************************************************************************************
		//The pair creation method works as follows.
		//we make a double loop over all particles i,j. 
		//for each particle we calculate the distance in between the two given particle using the nint function
		//the last equation helps us to keep the perioded boundary conditons in the simulation.
		//*****************************************************************************************************
		void updateFriends(){
			this->clearSimulationFriends();
			#pragma omp parallel for
			for(unsigned int i = 0; i < _amountobjects; ++i){
				#pragma omp parallel for
				for(unsigned int j = 0; j < i; ++j){
					if(i != j){	
						double xtemp = calcDistance(_objects[i].getPosition().getX(),_objects[j].getPosition().getX(),_systemsize);
						double ytemp = calcDistance(_objects[i].getPosition().getY(),_objects[j].getPosition().getY(),_systemsize);
						double ztemp = calcDistance(_objects[i].getPosition().getZ(),_objects[j].getPosition().getZ(),_systemsize);
						//now check wether the distance in between the given objects is smaller than rcutoff
						if(sqrt(xtemp*xtemp + ytemp*ytemp + ztemp*ztemp) < _rcutof){
							#pragma omp critical
							{
								_objects[i].addFriend(&_objects[j]);
								_objects[j].addFriend(&_objects[i]);
							}	
						}
					}
					
				}
			}
		}
	
	public:
		//*****************************************************************************
		//This function calculates <v(0)v(current)> being the velocity autocorrelation function.
		//*****************************************************************************
		double calcAutoCorrelation(){
			double currentguess = 0;
			vector<Vec> currentspeeds = this->getSpeeds();
			//we loop over all the particles in the system
			for(unsigned int i = 0; i<_amountobjects; ++i){
				double extra = _initialspeeds[i]*currentspeeds[i];
				currentguess += extra;		
			}
			//now all we have to do is to output this current result bij the amount of contributions being the amount of objects
			return currentguess/_amountobjects;
		}

	public:
		//***************************************************************
		//This function will loop over all the particles, calc the distance they have to neigbours and bin this data in the ouput
		//***************************************************************
		vector<double> binDistances(){
			//initialise an empty output vector
			vector<double> output(_subs,0.);
			#pragma omp parallel for
			for(unsigned int i = 0; i<_amountobjects; ++i){
				#pragma omp parallel for
				for(unsigned int j = 0; j<i; ++j){
					if(i!=j){
						//calculate the x,y,z distances between the objects
						double xtemp = calcDistance(_objects[i].getPosition().getX(),_objects[j].getPosition().getX(),_systemsize);
						double ytemp = calcDistance(_objects[i].getPosition().getY(),_objects[j].getPosition().getY(),_systemsize);
						double ztemp = calcDistance(_objects[i].getPosition().getZ(),_objects[j].getPosition().getZ(),_systemsize);
						//use this to calculate the distance in between these objects
						double distance = sqrt(xtemp*xtemp+ytemp*ytemp+ztemp*ztemp);
						//now, we check in which bin it belongs
						//eg. dist = 6.3 binsize = 3 --> floor(6.3/3)=floor(2.1)=bin 2.
						#pragma omp atomic
						++output[floor(distance/_subslength)];
					}
			}
		}
		return output;
	}
	
	public:
		//**********************************************************************
		//this function will calculate the total kinetic energy in the system.
		//**********************************************************************	
		double calcKineticEnergy(){
			double output = 0;
			for(unsigned int i = 0; i < _amountobjects; ++i){
				output += 0.5*pow(_objects[i].getSpeed().norm(),2);
			}
			return output;
		}
	public:
		//*********************************************************************
		//This function calculates the average speed in the system as in (sum_objets vx+vy+vz)/3*amountobjects
		//*********************************************************************
		double calcAverageSpeed(){
			double output = 0;
			for(unsigned int i = 0; i < _amountobjects; ++i){
				Vec speed = _objects[i].getSpeed();
				output += abs(speed.getX()) + abs(speed.getY()) + abs(speed.getZ());
			}
			return output/(3.*_amountobjects);
		}


	public:
		//************************************************************
		//This function will calculate the total potential energy in the system
		//***********************************************************
		double calcPotentialEnergy(){
			double output = 0;
			//first we get the energy due to all friend/friend interations. remember to devide the result by two !
			#pragma omp parallel for
			for(unsigned int i = 0; i< _amountobjects; ++i){
				Object on = _objects[i];
				#pragma omp parallel for
				for(unsigned int j = 0; j<on.getAmountFriends(); ++j ){
					Object by = on.getFriendAt(j);
					//*******
					//calculate the x,y,z distances between the two objects using the NINT function.
					//*******
					double xtemp = calcDistance(on.getPosition().getX(),by.getPosition().getX(),_systemsize);
					double ytemp = calcDistance(on.getPosition().getY(),by.getPosition().getY(),_systemsize);
					double ztemp = calcDistance(on.getPosition().getZ(),by.getPosition().getZ(),_systemsize);
					//****
					//calculate the actual r distance in between the objects
					//****
					double distance = sqrt(xtemp*xtemp+ytemp*ytemp+ztemp*ztemp);
					//calculate the potential between these two objects
					#pragma omp atomic
					output += 4*(pow(distance,-12)-pow(distance,-6));
				}
			}
			//at this point output contains double the potential of the friend friend interaction. we still need to take the non friend non friend into acount
			//this is done using the formule in the book 8.18
			//we get:
			output = output/2. + 8.*M_PI*(_amountobjects)*(_amountobjects)*(_amountobjects-1.)/(pow(_systemsize,3.))*(1./9.*pow(_rcutof,-9.)-1./3.*pow(_rcutof,-3.));
			return output;	
		}
	public:
		//*****************************************************************
		//This function will calculate a factor lambda(temp) so that lambda*v's behave according to the
		//equipartition theorem at that temperature
		//*****************************************************************
		double calcLambda(double temp){
			//the used formula is basicly sqrt(total energy/kinetic energy) 
			//the 119.8 is the depth of the lennard jones/boltzmann constant
			return sqrt((1.5*(_amountobjects-1)*temp/119.8)/(calcKineticEnergy()));
		}

	public:
		//**********************************************************************
		//This function will calculate the force on a given object
		//***********************************************************************
		vector<Vec> calcForce(){
			vector<Vec> output(_amountobjects,Vec());
			//we fire up a thread for every object in the simulation. Every thread calculates the force on that object and stores it
			#pragma omp parallel for
			for(unsigned int i=0; i<_amountobjects; ++i){
				//get the object under avaluation
				Object on = _objects[i];
				unsigned int amountfriends = on.getAmountFriends();
				//intialize the force on this object
				Vec forceon = Vec();
				#pragma omp parallel for
				for(unsigned int j = 0; j<amountfriends; ++j){
					//get the object that exterts the force on on.
					Object by = on.getFriendAt(j);
					//*******
					//calculate the x,y,z distances between the two objects using the NINT function.
					//*******
					double xtemp = calcDistance(on.getPosition().getX(),by.getPosition().getX(),_systemsize);
					double ytemp = calcDistance(on.getPosition().getY(),by.getPosition().getY(),_systemsize);
					double ztemp = calcDistance(on.getPosition().getZ(),by.getPosition().getZ(),_systemsize);
					//****
					//calculate the actual r distance in between the objects
					//****
					double distance = sqrt(xtemp*xtemp+ytemp*ytemp+ztemp*ztemp);
					//*******************
					//the force from i on j is calculated using
					//Fij_x = (xi-xj)*(48*rij^-14 - 24*rij^-8)
					//first we calculate the common postfactor for the x,y,z part
					//*******************
					double postfactor = 48*pow(distance,-14)-24*pow(distance,-8);
					//now we calculate the x,y and z parts of the force. We will store these in xtemp, ytemp,...
					xtemp = xtemp*postfactor;
					ytemp = ytemp*postfactor;
					ztemp = ztemp*postfactor;
					//add the new Vec to the outputi
					#pragma omp critical
					{
					forceon = forceon + Vec(xtemp,ytemp,ztemp);
					}
				}
				output[i] = forceon;
			}
			return output;
		}

};

//************************
//OPERERATOR OVERLOADING//
//***********************

//overloaded <<
//=============
ostream& operator <<(ostream& os,const Simulation& a){
	for(unsigned int i=0; i<a.getAmountobjects(); ++i){
		os << a.getTime() << '\t' << i << '\t' << a.getObjectAt(i) << endl;
	}
	return os;
}

//************************************************
//Function that makes a couple of verlet steps
//***********************************************
void simulate(Simulation& sim,double inittime,double time, double twriter, unsigned int scalesteps, double timestep, double temp){
	//get the gpu and walltime
	auto cpu0 = clock();
	auto wall0 = chrono::system_clock::now();
	ofstream coordinates("Data/coordinates.md");
	ofstream cond("Data/conditions.md");
	ofstream velcorr("Data/velocitycorr.md");	
	ofstream energy("Data/energy.md");	
	ofstream corr("Data/correlation.md");
	//the timestep inputted by the user is only a suggestion (of the correct order) we make a better guess using the density and the temperature
	//larger density --> smaller timestep, analog for temperature.
	//density is a r^3 effeect so pow(-1/3) while temp is related to speed which is a v^2 <-> r^2 effect.
	//taking these(rudementary) effects into consideration we get:
	timestep = timestep*pow(sim.getDensity(),-1./3.)*pow(2*temp,-1./2.);
	//using this timestepd we can calculate the amount of initsteps and simulatedsteps to be made
	unsigned int initsteps =(int)ceil(inittime/timestep);
	unsigned int simulatedsteps = (int)ceil(time/timestep);
	unsigned int writersteps = (int)ceil(twriter/timestep);

	cond << "#THe timestep of the simulation was: " << timestep << endl;
	cond << "#the amount of initialising time is: " << inittime << "*10^-4s"<< endl;
	cond << "#the actual simulated time is: " << time << "*10^-14s" << endl;
	cond << "#Every " << writersteps*timestep << "*10^-14s the data will be written away." << endl;
	cond << "#The amount of datasets in the file is: " << simulatedsteps/writersteps << endl;
	cond << "#the systemsize is: " << sim.getSystemsize() << " in units of sigma" << endl;
	cond << "#The temperature of the simulation was: " << temp << " K"  << endl;
	cond << "#THe density of the simulation was: " << sim.getDensity() << " in units of sigma^-3"  << endl;
	cond << "#There were " << sim.getAmountobjects() << " objects in the simulation" << endl;
	
	unsigned int currstep = 0;
	//we hold a vector in which we store the summed binned distances of the system.
	//in the end of the simulation this will give us <n(r)> which enables us to calculate the correlation function
	vector<double> summedbinned(sim.getSubs(),0.);
	unsigned int contributions = 0;
	//we hold a flag wether the initialspeeds have already been stored. Also, we store v^2 ini
	bool initialspeedsmade = false;
	double vsquaredini = 0;
	//we will calculate the next step at which to calculate the frienlist on as following:
	//<v(t)>*timestep*steps = rcutofextra where rcutofextra is a "safety zone" 
	//for now, the friendlist has to be updated at the first step!
	unsigned int updatefriendsstep = 0;
		
	//make the needed simulation steps
	cout << "Making the actual simulation" << endl;
	cout << "-------------------------------------------------------------------" << endl;
	while(currstep <= simulatedsteps+initsteps){
		cout << "Simulation progres: " << currstep <<"/"<<simulatedsteps+initsteps << '\r';
		cout.flush();
		//*************************
		//Init initialspeeds if needed
		//************************
		if(!initialspeedsmade && currstep >= initsteps){
			sim.set_Initialspeeds();
			vsquaredini = sim.calcKineticEnergy()*2/sim.getAmountobjects();
			initialspeedsmade = true;
		}

		//*******************
		//Update the friend list. Note that this is always done for currstep == 0
		//also, keep in mind that this is the only N^2 proces in the simulation. We try to avoid it as much as possible
		//******************
		if(currstep == updatefriendsstep){
			sim.updateFriends();
			int steps = (int)ceil(sim.getRcutofextra()/(timestep*5.*sim.calcAverageSpeed()));
				if(steps > 500){
					steps = 500;
					cout << "very large friendsteps" << endl;
				}
			updatefriendsstep += steps;
				
		}
		//******************************
		//rescale velocities when we are in init phase
		//***************************
		if(currstep%scalesteps == 0 && currstep < initsteps-writersteps){
			sim.multiplySpeeds(sim.calcLambda(temp));
		}
		//******************
		//Print away the system if needed
		//we start writing away before initsteps has ended so when the plot starts at t=0 it will look nicer...
		//in theory this value might be off a tiny bit as it is still in the init phase, but the error will be neglectable
		//*****************
		if(currstep%writersteps == 0 && currstep >= initsteps-writersteps){
			#pragma omp parallel num_threads(3)
			{
				if(omp_get_thread_num() == 0 ){
					//note that we insert two endlines.This is because gnuplot than sees the 2 blocks as different datasets!
					coordinates << sim << endl << endl;
				}
				if(omp_get_thread_num() == 1){	
					
					//we also want to output energy over time. This is done without the several blocks
					double potential = sim.calcPotentialEnergy();
					double kinetic = sim.calcKineticEnergy();
					energy << sim.getTime()-inittime << '\t' << potential << '\t' << kinetic << '\t' << potential + kinetic << endl;
				}
				if(omp_get_thread_num() == 2){
					//also we want to add this setup to the summedbinned vector
					summedbinned = summedbinned + sim.binDistances();
					contributions ++;
					//write away the velocity autocorrelation, note the normalisation
					velcorr << sim.getTime()-inittime << '\t' << sim.calcAutoCorrelation()/vsquaredini << endl;
				}
			}
		}	
		//********************
		//Make the actual step
		//********************
		//first you calculate the force array at the current situation
		vector<Vec> forcenow = sim.calcForce();
		//shift the positions of the particles using the just calculated array
		sim.shiftPositions(timestep*sim.getSpeeds()+timestep*timestep*0.5*forcenow);
		//calculate the force at this new position
		vector<Vec> forcesec = sim.calcForce();
		//shift the speeds
		sim.shiftSpeeds(timestep*0.5*(forcenow+forcesec));
		//update the current step
		currstep++;
		sim.shiftTime(timestep);
	}
	cout << endl << endl;
	//calculate <n(r)>
	summedbinned = (1./contributions)*summedbinned;
	//now we need to loop over all the elements in the summedbinned vector and multiply them with the correct prefactor
	//having done this we obtained the correclation function
	double prefactor = 2.*pow(sim.getSystemsize(),3.)/(sim.getAmountobjects()*(sim.getAmountobjects()-1.)*4.*M_PI);
	summedbinned = prefactor*summedbinned;
	//now every element gets an individual prefactor
	//after this is done the elements in the vector are the actual correlation functions so we write them outan@Gertian-Pc:~/CourseNotes/Computationele/Code/MD/Data$ ffmpeg -framerate 25 -i Plot%04d.png -c:v libx264 -r 30 -pix_fmt bgr565 out.mp4

	unsigned int subs = sim.getSubs();
	for(unsigned int i = 0; i<subs; ++i){
		prefactor = 1./(pow((i+1.)*sim.getSubslength(),2.)*sim.getSubslength());
		corr << (i+1.)*sim.getSubslength() << "\t" <<summedbinned[i]*prefactor << endl; 
	}

	auto cpu1 = clock();
	auto wall1 = chrono::system_clock::now();

	//NOW some rather ugly non windows compatible code.
	//For non linux users, do the following, 
	//1) run the code and see if it works, it might I don't know...
	//2) if it does not work, comment this codeblock(all the system... commands) out and recompile
	//3) Check in the map Data, there is a file 'conditions.md' in this file you will find the systemsize and the amount of datapoints made
	//4) in the main directory, run gnuplot -e "sytemwidth=A; plots=B" Plots.gnu where A and B are the systemsize and the amount of datapoints
	//5) Finally check in the map Data, the plotted data can be found there
	cout.flush();
	cout << "ploting all the generated data" << endl; 
	cout << "-------------------------------------------------------------------" << endl;
	string order = "gnuplot -e \"systemwidth="+to_string(sim.getSystemsize())+";plots="+to_string(simulatedsteps/writersteps-1.)+";highlight="+to_string(ceil(sim.getAmountobjects()/2.))+";tfinal="+to_string(time)+"\" Plots.gnu";
	system(order.c_str());
	cout << endl << endl;
	cout << "Making a call to ffmpeg to make a movie out of the simulation" << endl;
	cout << "-------------------------------------------------------------------" << endl;
	system("ffmpeg -y -framerate 20 -i ./Data/Plot%04d.png -c:v libx264 -r 30 -pix_fmt bgr565 ./Data/out.mp4");
	cout << endl << endl;
	system("rm ./Data/Plot*.png");
	cout << "done" << endl;
	cout << "-------------------------------------------------------------------" << endl;
	cout << "used CPU  time is: " << (double)(cpu1 - cpu0)/CLOCKS_PER_SEC << "s" << endl;	
	cout << "used wall time is: " << (double)(wall1 - wall0).count()*(1./1000000000.) << "s"  << endl;
}

//***********************************************
//***********************************************
//
//
//The main loop of the program
//
//
//***********************************************
//***********************************************
int main(){
	//make a simulation with 1 as density. THis means 1 particle per volume in units of potential length 
	//3*3 = 9 boxes in the simulation and 3*sigma as the rcutof.
	//Also, we will make 100 subdivisions for the calculation of the correlation length
	cout << "Gathering input: " << endl;
	cout << "-------------------------------------------------------------------" << endl;
	double density;
	cout << "Please enter the density to run the simulation at" << endl;
	cin >> density;
	
	int cubes;
	cout << "How many equivalent systems do you want to include in the system ?" << endl;
	cin >> cubes;
	
	double inittime;
	cout << "How long should the init run ? Keep in mind: Order(timestep = 0.004)" << endl;
	cin >> inittime;

	double time;
	cout << "How long should the simulation run ? Keep in mind: Order(timestep = 0.004)" << endl;
	cin >> time; 

	double temp;
	cout << "please enter the temperature to run the simulation at" << endl;
	cin >> temp;
	
	double twriter;
	cout << "please enter the time you want to have in between 'datapoints'" << endl;
	cin >> twriter;
	cout << endl;
	cout << endl;

	Simulation a = Simulation(density,cubes,2.5,0.4,320);
	simulate(a,inittime,time,twriter,5,0.001,temp);

	return 0;
}
