/*program to simulate the two dimensional ising model in equilibrium with a heat bath; Kb=1.0*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
void main()
{
long unsigned seed;
long nspin, inspin, jnspin, spin[100][100], accept, nconfig, in_config, s_accept, knspin, lnspin, spin_, mcs, imcs;
double rnd, rnd1, rnd2, rnd3, energy_old, energy_trial, sigma1, sigma2, sigma3, in_j, std_dev_, free_e, entropy;
double avg_e, avg2_e, energy_accum, energy2_accum, del_e, beta, temp, var, energy, prob, probab, wp;
FILE *infile;
FILE *outfile;
infile=fopen("ising2d.in","r");
outfile=fopen("ising2d.out","w");
fscanf(infile,"%ld",&spin_);
fprintf(outfile,"Number of particles in each dimension [INPUT]:%ld \n", spin_);
nspin=spin_*spin_;
fprintf(outfile,"Number of particles in 2D [INPUT]:%ld\n",nspin);
fscanf(infile,"%ld",&mcs);
fscanf(infile,"%ld",&nconfig);
fscanf(infile,"%lf",&in_j);
fprintf(outfile,"Interaction coefficient [INPUT]:%lf\n",in_j);
fscanf(infile,"%lf",&temp);
fprintf(outfile,"Temperature [INPUT]:%lf\n",temp);
beta=1.0/temp;
seed=time(NULL);
srand(seed);
for(inspin=0;inspin<spin_;inspin++) //Initialize configuration
for(jnspin=0;jnspin<spin_;jnspin++)
{
rnd1= rand()/(RAND_MAX+1.0);
if(rnd1>0.5)
spin[inspin][jnspin]=1;
else
spin[inspin][jnspin]=-1;
}
sigma1=0.0;
for(inspin=0;inspin<spin_;inspin++) //Energy of initial configuration
for(jnspin=0;jnspin<spin_;jnspin++)
{
spin[spin_][0]=spin[0][0];
spin[0][spin_]=spin[0][0];
spin[inspin][spin_]=spin[inspin][0];
spin[spin_][jnspin]=spin[0][jnspin];
sigma1= sigma1+ (spin[inspin][jnspin] *spin[inspin+1][jnspin]) +(spin[inspin][jnspin]*spin[inspin][jnspin+1]);
}
energy_old=-in_j*(sigma1/2.0);
energy2_accum=0.0;
energy_accum=0.0;
accept=0;
for(inspin=0;inspin<spin_;inspin++) //main loop
for(jnspin=0;jnspin<spin_;jnspin++)
{
if(spin[inspin][jnspin]==1) //make trial change in spin
spin[inspin][jnspin]=(-1);
else
spin[inspin][jnspin]=1;
sigma2=0.0;
for(knspin=0;knspin<spin_;knspin++) //energy after trial change
for(lnspin=0;lnspin<spin_;lnspin++)
{
spin[spin_][0]=spin[0][0];
spin[0][spin_]=spin[0][0];
spin[knspin][spin_-1]=spin[knspin][0];
spin[spin_-1][jnspin]=spin[0][jnspin];
sigma2=sigma2+(spin[knspin][lnspin]* spin[knspin+1] [lnspin])+(spin[knspin][lnspin]*spin[knspin][lnspin+1]);
}
energy_trial=-in_j*(sigma2/2.0);
del_e=energy_trial-energy_old; //finding change in energy due random trial
if(del_e<=0.0) //accept the change
{
energy_accum=energy_accum+energy_trial;
energy2_accum+=pow(energy_trial,2.0);
accept++;
energy_old=energy_trial;
}
else
{
rnd3 = rand()/(RAND_MAX+1.0);
if(rnd3<exp(-beta*del_e)) //accept the change
{
energy_accum=energy_accum+energy_trial;
energy2_accum+=pow(energy_trial,2.0);
accept++;
energy_old=energy_trial;
}
else
{
energy_old=energy_old; //change is not accepted
if(spin[inspin][jnspin]==1) //changing the configuration to the initial state
spin[inspin][jnspin]=-1;
else
spin[inspin][jnspin]=1;
}
}
}
fprintf(outfile,"Total number of configurations= %ld\n", spin_*spin_);
fprintf(outfile,"Accepted configurations for energy calculation =%ld\n",accept); //number of accepted moves
avg_e=energy_accum/(double)accept; //mean energy of the system
avg2_e=energy2_accum/(double)accept;
fprintf(outfile,"Mean Energy=%lf(K)\n",avg_e);
var=avg2_e-pow(avg_e,2.0); //variance
fprintf(outfile,"Standard deviation of energy values=%lf \n\n", sqrt(var));
s_accept=0;
for(imcs=0;imcs<mcs;imcs++) //entropy calculation
{
for(in_config=0;in_config<nconfig;in_config++) //number of configuration
{
for(inspin=0;inspin<spin_;inspin++) //assigning random spin
for(jnspin=0;jnspin<spin_;jnspin++)
{
rnd= rand()/(RAND_MAX+1.0);
if(rnd>0.5)
spin[inspin][jnspin]=1;
else
spin[inspin][jnspin]=-1;
}
sigma3=0.0;
for(knspin=0;knspin<spin_;knspin++) //Energy calculation
for(lnspin=0;lnspin<spin_;lnspin++)
{
spin[spin_][0]=spin[0][0];
spin[0][spin_]=spin[0][0];
spin[spin_-1][lnspin]=spin[0][lnspin];
spin[knspin][spin_-1]=spin[knspin][0];
sigma3=sigma3+(spin[knspin][lnspin]* spin[knspin+1]
[lnspin])+(spin[knspin][lnspin]*spin[knspin][lnspin+1]);
}
energy=-in_j*(sigma3/2.0);
if((energy>(avg_e-sqrt(var)))&&(energy<(avg_e+sqrt(var)))) //accept the microstate
s_accept++;
}
prob=(double)s_accept/(double)nconfig;
probab+=prob;
}
wp=(double)probab/(double)(mcs*s_accept);
fprintf(outfile,"Number of Monte Carlo steps [INPUT]:%ld\n",mcs);
fprintf(outfile,"Number of initial configurations for entropy [INPUT]:%ld\n",nconfig);
fprintf(outfile,"Total number of configurations= %ld\n", mcs*nconfig);
fprintf(outfile,"Accepted configurations for entropy calculation=%ld\n",s_accept); //number of accepted moves
fprintf(outfile,"Probability=%lf\n",wp); //acceptance probability
entropy=log(wp); //s=kb ln W
fprintf(outfile,"Entropy=%lf(K)\n",entropy);
free_e=avg_e-temp*entropy; //g=e-t*s
fprintf(outfile,"Free energy=%lf(K)\n",free_e);
fclose(infile);
fclose(outfile);
}