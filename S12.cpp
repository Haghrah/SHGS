#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>
//===================================================================================
using namespace std;
//===================================================================================
#define hydro_num 4
#define thermo_num 1
#define population_num 500
#define mutation_num 2000
#define GNG 1000
//C at the end of definitions stands for constant
//L & H ................................ Low & High
#define TCC 0.000001
#define PTLC 250
#define PTHC 250
#define QHLC 250
#define QHHC 250
#define VHLC 250
#define VHHC 250
#define VENDC 250
#define PHLC 250
#define PHHC 250
#define VINITC 250

//===================================================================================
//HYDRO_DATA START
double C_hydro[hydro_num][6]={-0.0042,-0.42,0.03,0.9,10,-50,
							-0.004,-0.3,0.015,1.14,9.5,-70,
							-0.0016,-0.3,0.014,0.55,5.5,-40,
							-0.003,-0.31,0.027,1.44,14,-90};
double V_hydro_min[hydro_num]={80,60,100,70};
double V_hydro_max[hydro_num]={150,120,240,160};
double V_hydro_init[hydro_num]={100,80,170,120};
double V_hydro_end[hydro_num]={120,70,170,140};
double Q_hydro_min[hydro_num]={5,6,10,13};
double Q_hydro_max[hydro_num]={15,15,30,25};
double P_hydro_min[hydro_num]={0,0,0,0};
double P_hydro_max[hydro_num]={500,500,500,500};
int Time_delay_hydro[hydro_num]={2,3,4,0};
int K_hydro[hydro_num][hydro_num-1]={-1,-1,-1,
									-1,-1,-1,
									0,1,-1,
									2,-1,-1};
double inflow_t[hydro_num][24]={10,9,8,7,6,7,8,9,10,11,12,10,11,12,11,10,9,8,7,6,7,8,9,10,
								8,8,9,9,8,7,6,7,8,9,9,8,8,9,9,8,7,6,7,8,9,9,8,8,
								8.1,8.2,4,2,3,4,3,2,1,1,1,2,4,3,3,2,2,2,1,1,2,2,1,0,
								2.8,2.4,1.6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double S_t[hydro_num][24]={0};
//HYDRO_DATA END
//===================================================================================
//THERMO_DATA START
double C_thermo[thermo_num][5]={0.002,19.2,5000,700,0.085};
double P_thermo_min[thermo_num]={500};
double P_thermo_max[thermo_num]={2500};
//THERMO_DATA END
//===================================================================================
double P_load[24]={1370,1390,1360,1290,1290,1410,1650,2000,2240,2320,2230,2310,2230,2200,2130,2070,2130,2140,2240,2280,2240,2120,1850,1590
};
//===================================================================================
struct point{
	double Q_hydro_negative[hydro_num][24],Q_hydro[hydro_num][24],V_hydro[hydro_num][24],P_hydro[hydro_num][24],P_thermo[thermo_num][24];
	double fitness;
}population[population_num];
//===================================================================================
double Q_hydro(point in_point,int number,int time);
//===================================================================================
int sgn(double d){ 
    return d<0?-1:d>0;
}
//===================================================================================
int compare (const void * in1, const void * in2)
{
	if(((point*)in1)->fitness<((point*)in2)->fitness) return -1;
	else return 1;
}
//===================================================================================
double S_hydro(point in_point,int number,int time)
{
	return 0;
}
//===================================================================================
double V_hydro(point in_point,int number,int time)
{
	if(time==-1){
		return V_hydro_init[number];
	}else{
		return in_point.V_hydro[number][time];
	}
}
//===================================================================================
double Q_hydro(point in_point,int number,int time)
{
	if(time<0)
	{
		return in_point.Q_hydro_negative[number][24+time];
	}else{
		double output=0;
		int i;
		output+=(V_hydro(in_point,number,time-1)+inflow_t[number][time]-in_point.V_hydro[number][time]-S_t[number][time]);
		for(i=0;i<hydro_num-1;i++)
		{
			if(K_hydro[number][i]!=-1)
			{
				output+=Q_hydro(in_point,K_hydro[number][i],time-Time_delay_hydro[K_hydro[number][i]]);
				output+=S_hydro(in_point,K_hydro[number][i],time-Time_delay_hydro[K_hydro[number][i]]);
			}
		}
		return output;
	}
}
//===================================================================================
double P_hydro(point in_point,int number,int time)
{
	return (C_hydro[number][0]*in_point.V_hydro[number][time]*in_point.V_hydro[number][time]+C_hydro[number][1]*in_point.Q_hydro[number][time]*in_point.Q_hydro[number][time]+C_hydro[number][2]*in_point.V_hydro[number][time]*in_point.Q_hydro[number][time]+C_hydro[number][3]*in_point.V_hydro[number][time]+C_hydro[number][4]*in_point.Q_hydro[number][time]+C_hydro[number][5]);
}
//===================================================================================
double P_load_point(point in_point,int time)
{
	double sum=0;
	int i;
	for(i=0;i<hydro_num;i++)
		sum+=in_point.P_hydro[i][time];
	for(i=0;i<thermo_num;i++)
		sum+=in_point.P_thermo[i][time];

	return sum;
}
//===================================================================================
double TC(point in_point)
{
	int i,j;
	double output=0;
	for(i=0;i<thermo_num;i++)
	{
		for(j=0;j<24;j++)
		{
			output+=(C_thermo[i][0]*in_point.P_thermo[i][j]*in_point.P_thermo[i][j]+C_thermo[i][1]*in_point.P_thermo[i][j]+C_thermo[i][2]+abs(C_thermo[i][3]*sin(C_thermo[i][4]*(P_thermo_min[i]-in_point.P_thermo[i][j]))));
		}
	}
	return output;
}
//===================================================================================
double fitness(point in_point)
{
	int i,j;
	double output=0;
	//TC
	output+=TC(in_point)*TCC;
	//P_thermo
	for(i=0;i<thermo_num;i++)
	{
		for(j=0;j<24;j++)
		{
			if(in_point.P_thermo[i][j]<P_thermo_min[i]) 
			{
				output+=(pow(abs(in_point.P_thermo[i][j]-P_thermo_min[i]),1)+PTLC*abs(in_point.P_thermo[i][j]-P_thermo_min[i]));
			}else if(in_point.P_thermo[i][j]>P_thermo_max[i])
			{
				output+=(pow(abs(in_point.P_thermo[i][j]-P_thermo_max[i]),1)+PTHC*abs(in_point.P_thermo[i][j]-P_thermo_max[i]));
			}
		}
	}
	//Q_hydro
	for(i=0;i<hydro_num;i++)
	{
		for(j=0;j<24;j++)
		{
			if(in_point.Q_hydro[i][j]<Q_hydro_min[i])
			{
				output+=(pow(abs(in_point.Q_hydro[i][j]-Q_hydro_min[i]),5)+QHLC*abs(in_point.Q_hydro[i][j]-Q_hydro_min[i]));
			}else if(in_point.Q_hydro[i][j]>Q_hydro_max[i])
			{
				output+=(pow(abs(in_point.Q_hydro[i][j]-Q_hydro_max[i]),5)+QHHC*abs(in_point.Q_hydro[i][j]-Q_hydro_max[i]));
			}
		}
	}

	//Q_hydro_negative
	for(i=0;i<hydro_num;i++)
	{
		for(j=0;j<24;j++)
		{
			if(in_point.Q_hydro_negative[i][j]<Q_hydro_min[i])
			{
				output+=(pow(abs(in_point.Q_hydro_negative[i][j]-Q_hydro_min[i]),5)+QHLC*abs(in_point.Q_hydro_negative[i][j]-Q_hydro_min[i]));
			}else if(in_point.Q_hydro_negative[i][j]>Q_hydro_max[i])
			{
				output+=(pow(abs(in_point.Q_hydro_negative[i][j]-Q_hydro_max[i]),5)+QHHC*abs(in_point.Q_hydro_negative[i][j]-Q_hydro_max[i]));
			}
		}
	}
	
	//V_hydro
	for(i=0;i<hydro_num;i++)
	{
		for(j=0;j<24;j++)
		{
			if(in_point.V_hydro[i][j]<V_hydro_min[i])
			{
				output+=(pow(abs(in_point.V_hydro[i][j]-V_hydro_min[i]),1)+VHLC*abs(in_point.V_hydro[i][j]-V_hydro_min[i]));
			}else if(V_hydro(in_point,i,j)>V_hydro_max[i])
			{
				output+=(pow(abs(in_point.V_hydro[i][j]-V_hydro_max[i]),1)+VHHC*abs(in_point.V_hydro[i][j]-V_hydro_max[i]));
			}
		}
	}


	//V[24]=V_End
	for(i=0;i<hydro_num;i++)
	{
		output+=(pow(abs(in_point.V_hydro[i][23]-V_hydro_end[i]),1)+VENDC*abs(in_point.V_hydro[i][23]-V_hydro_end[i]));
	}

	//V[0]=V_init
	for(i=0;i<hydro_num;i++)
	{
		output+=(pow(abs(V_hydro(in_point,i,-1)-V_hydro_init[i]),1)+VINITC*abs(V_hydro(in_point,i,-1)-V_hydro_init[i]));
	}


	//P_hydro
	for(i=0;i<hydro_num;i++)
	{
		for(j=0;j<24;j++)
		{
			if(in_point.P_hydro[i][j]<P_hydro_min[i])
			{
				output+=(pow(abs(in_point.P_hydro[i][j]-P_hydro_min[i]),1)+PHLC*abs(in_point.P_hydro[i][j]-P_hydro_min[i]));
			}else if(P_hydro(in_point,i,j)>P_hydro_max[i])
			{
				output+=(pow(abs(in_point.P_hydro[i][j]-P_hydro_max[i]),1)+PHHC*abs(in_point.P_hydro[i][j]-P_hydro_max[i]));
			}
		}
	}
	
	return output;
}
//===================================================================================
point generate_point()
{
	point out_point;
	int i,j,k;
	for(i=0;i<hydro_num;i++)
	{
		for(j=0;j<24;j++)
		{
			out_point.V_hydro[i][j]=double(rand())/RAND_MAX*(V_hydro_max[i]-V_hydro_min[i])+V_hydro_min[i];
			out_point.Q_hydro_negative[i][j]=double(rand())/RAND_MAX*(Q_hydro_max[i]-Q_hydro_min[i])+Q_hydro_min[i];
		}
		out_point.V_hydro[i][23]=V_hydro_end[i];
	}
	for(i=0;i<hydro_num;i++)
	{
		for(j=0;j<24;j++)
		{
			out_point.Q_hydro[i][j]=Q_hydro(out_point,i,j);
			out_point.P_hydro[i][j]=P_hydro(out_point,i,j);
		}
	}
	for(i=0;i<thermo_num;i++)
	{
		for(j=0;j<24;j++)
		{
			out_point.P_thermo[i][j]=P_load[j];
			for(k=0;k<hydro_num;k++)
			{
				out_point.P_thermo[i][j]-=out_point.P_hydro[k][j];
			}
		}
	}
	out_point.fitness=fitness(out_point);
	return out_point;
}
//===================================================================================
void generate_initial_population(point *points,int point_count)
{
	int i;
	for(i=0;i<point_count;i++)
	{
		points[i]=generate_point();
	}
	qsort(points,population_num,sizeof(point),compare);
}
//===================================================================================
point combination(point chromo1,point chromo2)
{
	point outpoint;
	int i,j,k;
	double c1,c2;
	for(i=0;i<hydro_num;i++)
	{
		for(j=0;j<24;j++)
		{
			c1=25*double(rand())/RAND_MAX;
			c2=5*double(rand())/RAND_MAX;
			outpoint.V_hydro[i][j]=(chromo1.V_hydro[i][j]*c1+chromo2.V_hydro[i][j]*c2)/(c1+c2);
			outpoint.Q_hydro_negative[i][j]=(chromo1.Q_hydro_negative[i][j]*c1+chromo2.Q_hydro_negative[i][j]*c2)/(c1+c2);
		}
		outpoint.V_hydro[i][23]=V_hydro_end[i];
	}
	for(i=0;i<hydro_num;i++)
	{
		for(j=0;j<24;j++)
		{
			outpoint.Q_hydro[i][j]=Q_hydro(outpoint,i,j);
			outpoint.P_hydro[i][j]=P_hydro(outpoint,i,j);
		}
	}
	for(i=0;i<thermo_num;i++)
	{
		for(j=0;j<24;j++)
		{
			outpoint.P_thermo[i][j]=P_load[j];
			for(k=0;k<hydro_num;k++)
			{
				outpoint.P_thermo[i][j]-=outpoint.P_hydro[k][j];
			}
		}
	}
	outpoint.fitness=fitness(outpoint);
	return outpoint;
}

long double gama(int CurrentGeneration)
{
	int i,a,b;
	a=8*double(CurrentGeneration)/GNG;
	b=a+4;
	long double output=0;
	for(i=a;i<b;i++)
		output+=0.5*(1+sgn(double(rand())/RAND_MAX-1+(i-a)/8.0))*pow((long double)(2),-i);
	return output;
}
point mutate(point in_point,int CurrentGeneration)
{
	point outpoint;
	int i,j,k;
	for(i=0;i<hydro_num;i++)
	{
		for(j=0;j<24;j++)
		{
			outpoint.V_hydro[i][j]=in_point.V_hydro[i][j]+sgn(double(rand())/RAND_MAX-0.5)*0.1*(V_hydro_max-V_hydro_min)*gama(CurrentGeneration);
			outpoint.Q_hydro_negative[i][j]=in_point.Q_hydro_negative[i][j]+sgn(double(rand())/RAND_MAX-0.5)*0.1*(Q_hydro_max-Q_hydro_min)*gama(CurrentGeneration);
		}
		outpoint.V_hydro[i][23]=V_hydro_end[i];
	}
	for(i=0;i<hydro_num;i++)
	{
		for(j=0;j<24;j++)
		{
			outpoint.Q_hydro[i][j]=Q_hydro(outpoint,i,j);
			outpoint.P_hydro[i][j]=P_hydro(outpoint,i,j);
		}
	}
	for(i=0;i<thermo_num;i++)
	{
		for(j=0;j<24;j++)
		{
			outpoint.P_thermo[i][j]=P_load[j];
			for(k=0;k<hydro_num;k++)
			{
				outpoint.P_thermo[i][j]-=outpoint.P_hydro[k][j];
			}
		}
	}
	outpoint.fitness=fitness(outpoint);
	return outpoint;
}

//===================================================================================
void generate_next_generation(point *points,int point_count,int CurrentGeneration)
{
	int i,j,k;
	point offspring;
	point mutated;
	for(i=0;i<(population_num/2);i++)
	{
		offspring=combination(points[i],points[population_num-i-1]);
		if(offspring.fitness<points[population_num-i-1].fitness)
		{
			points[population_num-i-1]=offspring;
		}
	}
	qsort(points,population_num,sizeof(point),compare);
	for(i=0;i<mutation_num;i++)
	{
		j=int(double(rand())/RAND_MAX*(population_num-1))/10;
		mutated=mutate(points[j],CurrentGeneration);
		if(mutated.fitness<points[population_num-1].fitness)
		{
			//points[population_num-1]=mutated;
			k=population_num-1;
			while( (points[k].fitness>mutated.fitness) && (k!=0) )
			{
				points[k]=points[k-1];
				k--;
			}
			points[k]=mutated;
			//qsort(points,population_num,sizeof(point),compare);
		}
	}
}

//===================================================================================
double TC_average(point *points,int point_count)
{
	int i;
	long double output=0;
	for(i=0;i<point_count;i++)
	{
		output+=TC(points[i]);
	}
	output/=point_count;
	return output;
}
//===================================================================================
int make_output(point in_point)
{
	ofstream outfile;
	outfile.open("outfile.html",ios_base::app);
	int i,j;
	if(!outfile)
	{
		cout<<"Cannot open output file.\n";
		return 1;
	}else{
		outfile.precision(5);
		outfile<<fixed;
		outfile<<"<html>\n<head>\n<title>Economic Power Dispatch Optimization Project</title>\n";
		outfile<<"</head>\n<body style=\"font-family:Courier;font-size:12px;\">\n";
		outfile<<"Hydro Data :\n<br/>\n<br/>\n<table border=\"1\" style=\"font-family:Courier;font-size:13px;\">\n";
		//-------------------------------------------------------------
		outfile<<"<tr>\n<th scope=\"col\">Time</th>\n";
		for(i=-1;i<24;i++)
		{
			outfile<<"<th scope=\"col\">"<<i+1<<"</th>\n";
		}
		outfile<<"</tr>\n";
		//-------------------------------------------------------------
		for(i=0;i<hydro_num;i++)
		{
			outfile<<"<tr>\n<th scope=\"col\">P_Hydro["<<i+1<<"]</th>\n";
			outfile<<"<th scope=\"col\">-</th>\n";
			for(j=0;j<24;j++)
			{
				outfile<<"<th scope=\"col\"";
				if(in_point.P_hydro[i][j]>P_hydro_max[i] || in_point.P_hydro[i][j]<P_hydro_min[i])
					outfile<<" bgcolor=\"#F87882\"";
				else
					outfile<<" bgcolor=\"#66FF66\"";
				outfile<<">"<<in_point.P_hydro[i][j]<<"</th>\n";
			}
			outfile<<"</tr>\n";
		}
		//-------------------------------------------------------------
		for(i=0;i<hydro_num;i++)
		{
			outfile<<"<tr>\n<th scope=\"col\">Q_Hydro["<<i+1<<"]</th>\n";
			outfile<<"<th scope=\"col\">-</th>\n";
			for(j=0;j<24;j++)
			{
				outfile<<"<th scope=\"col\"";
				if(in_point.Q_hydro[i][j]>Q_hydro_max[i] || in_point.Q_hydro[i][j]<Q_hydro_min[i])
					outfile<<" bgcolor=\"#F87882\"";
				else
					outfile<<" bgcolor=\"#66FF66\"";
				outfile<<">"<<in_point.Q_hydro[i][j]<<"</th>\n";
			}
			outfile<<"</tr>\n";
		}
		for(i=0;i<hydro_num;i++)
		{
			outfile<<"<tr>\n<th scope=\"col\">Q_Hydro_negative["<<i+1<<"]</th>\n";
			outfile<<"<th scope=\"col\">-</th>\n";
			for(j=0;j<24;j++)
			{
				outfile<<"<th scope=\"col\"";
				if(in_point.Q_hydro_negative[i][j]>Q_hydro_max[i] || in_point.Q_hydro_negative[i][j]<Q_hydro_min[i])
					outfile<<" bgcolor=\"#F87882\"";
				else
					outfile<<" bgcolor=\"#66FF66\"";
				outfile<<">"<<in_point.Q_hydro_negative[i][j]<<"</th>\n";
			}
			outfile<<"</tr>\n";
		}
		//-------------------------------------------------------------
		for(i=0;i<hydro_num;i++)
		{
			outfile<<"<tr>\n<th scope=\"col\">V_Hydro["<<i+1<<"]</th>\n";
			for(j=-1;j<24;j++)
			{
				outfile<<"<th scope=\"col\"";
				if(V_hydro(in_point,i,j)>V_hydro_max[i] || V_hydro(in_point,i,j)<V_hydro_min[i])
					outfile<<" bgcolor=\"#F87882\"";
				else
					outfile<<" bgcolor=\"#66FF66\"";
				outfile<<">"<<V_hydro(in_point,i,j)<<"</th>\n";
			}
			outfile<<"</tr>\n";
		}
		//-------------------------------------------------------------
		for(i=0;i<thermo_num;i++)
		{
			outfile<<"<tr>\n<th scope=\"col\">P_Thermo["<<i+1<<"]</th>\n";
			outfile<<"<th scope=\"col\">-</th>\n";
			for(j=0;j<24;j++)
			{
				outfile<<"<th scope=\"col\"";
				if(in_point.P_thermo[i][j]>P_thermo_max[i] || in_point.P_thermo[i][j]<P_thermo_min[i])
					outfile<<" bgcolor=\"#F87882\"";
				else
					outfile<<" bgcolor=\"#66FF66\"";
				outfile<<">"<<in_point.P_thermo[i][j]<<"</th>\n";
			}
			outfile<<"</tr>\n";
		}
		//-------------------------------------------------------------
		outfile<<"<tr>\n<th scope=\"col\">P_Load</th>\n";
		outfile<<"<th scope=\"col\">-</th>\n";
		for(j=0;j<24;j++)
		{
			outfile<<"<th scope=\"col\">"<<P_load_point(in_point,j)<<"</th>\n";
		}
		outfile<<"</tr>\n";
		//-------------------------------------------------------------
		outfile<<"</table>\n<br/><br/>";
		outfile<<"TC : "<<TC(in_point)<<"\n<br/><br/>";
		outfile<<"Fitness : "<<fitness(in_point)<<"\n<br/><br/>";
		outfile<<"<br/><br/><br/>\n</body>\n</html>";
	}
	outfile.close();
	return 0;
}
//===================================================================================
int make_raw_output(point in_point)
{
	ofstream outfile;
	outfile.open("input.txt",ios_base::app);
	int i,j;
	if(!outfile)
	{
		cout<<"Cannot open output file.\n";
		return 1;
	}
	outfile.precision(10);
	outfile<<fixed;
	outfile<<"P_hydro=[";
	for(i=0;i<hydro_num;i++)
	{
		outfile<<"[";
		for(j=0;j<24;j++)
		{
			outfile<<in_point.P_hydro[i][j]<<" ";
		}
		outfile<<"]\n";
	}
	outfile<<"];\n";

	outfile<<"Q_hydro=[";
	for(i=0;i<hydro_num;i++)
	{
		outfile<<"[";
		for(j=0;j<24;j++)
		{
			outfile<<in_point.Q_hydro[i][j]<<" ";
		}
		outfile<<"]\n";
	}
	outfile<<"];\n";

	outfile<<"Q_hydro_negative=[";
	for(i=0;i<hydro_num;i++)
	{
		outfile<<"[";
		for(j=0;j<24;j++)
		{
			outfile<<in_point.Q_hydro_negative[i][j]<<" ";
		}
		outfile<<"]\n";
	}
	outfile<<"];\n";

	outfile<<"V_hydro=[";
	for(i=0;i<hydro_num;i++)
	{
		outfile<<"[";
		for(j=0;j<24;j++)
		{
			outfile<<in_point.V_hydro[i][j]<<" ";
		}
		outfile<<"]\n";
	}
	outfile<<"];\n";

	outfile<<"P_thermo=[";
	for(i=0;i<thermo_num;i++)
	{
		outfile<<"[";
		for(j=0;j<24;j++)
		{
			outfile<<in_point.P_thermo[i][j]<<" ";
		}
		outfile<<"]\n";
	}
	outfile<<"];\n";

	outfile.close();
	return 0;
}
//===================================================================================
int make_TC_out(point *points,int point_count)
{
	ofstream outfile1;
	ofstream outfile2;
	ofstream outfile3;
	outfile1.open("TC_min.txt",ios_base::app);
	if(!outfile1)
	{
		cout<<"Cannot open output file.\n";
		return 1;
	}
	outfile1.precision(10);
	outfile1<<fixed;
	outfile1<<TC(points[0])<<" ";
	outfile1.close();

	outfile2.open("TC_max.txt",ios_base::app);
	if(!outfile2)
	{
		cout<<"Cannot open output file.\n";
		return 1;
	}
	outfile2.precision(10);
	outfile2<<fixed;
	outfile2<<TC(points[point_count-1])<<" ";
	outfile2.close();

	outfile3.open("TC_average.txt",ios_base::app);
	if(!outfile3)
	{
		cout<<"Cannot open output file.\n";
		return 1;
	}
	outfile3.precision(10);
	outfile3<<fixed;
	outfile3<<TC_average(points,point_count)<<" ";
	outfile3.close();
	return 0;
}
//===================================================================================
int make_TC_average_out(point *points,int point_count)
{
	ofstream outfile;
	outfile.open("TCAverage.txt",ios_base::app);
	if(!outfile)
	{
		cout<<"Cannot open output file.\n";
		return 1;
	}
	outfile.precision(10);
	outfile<<fixed;
	outfile<<TC_average(points,point_count)<<" ";
	outfile.close();
	return 0;
}
//===================================================================================
int main()
{
	clock_t t0=clock();
	int i/*,j*/;
	
	cout.precision(10);
	cout<<fixed;
	srand(static_cast<unsigned>(time(0)));
	generate_initial_population(population,population_num);
	
	for(i=1;i<=GNG;i++)
	{
		
		generate_next_generation(population,population_num,i);
		make_TC_out(population,population_num);
		make_output(population[0]);
		/*
		if(i%1000==0)
		{
			cin>>j;
			if(j=0)
				break;
		}
		*/
		cout<<"#"<<i<<"\n";
		cout<<"Fitness : "<<population[0].fitness<<"\n";
		cout<<"TC : "<<TC(population[0])<<"\n\n";
		
	}
	//make_TC_out(population,population_num);
	//make_output(population[0]);
	make_raw_output(population[0]);
	cout<<(float(clock()-t0))/CLOCKS_PER_SEC;
	cin>>i;
	return 0;
}
