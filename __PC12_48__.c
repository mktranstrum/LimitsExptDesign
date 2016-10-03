#include <math.h>
#include <stdio.h>
#include <float.h>
#define exponentiale M_E
#define pi M_PI
double max(double a, double b){
return a > b ? a : b;}
double min(double a, double b){
return a < b ? a : b;}
double root(double n,double x);

double root_0(double n,double x);

double root_1(double n,double x);

double cot(double x);

double cot_0(double x);

double arccot(double x);

double arccot_0(double x);

double coth(double x);

double coth_0(double x);

double csc(double x);

double csc_0(double x);

double arccsc(double x);

double arccsc_0(double x);

double csch(double x);

double csch_0(double x);

double sec(double x);

double sec_0(double x);

double arcsec(double x);

double arcsec_0(double x);

double sech(double x);

double sech_0(double x);

void var_types( int *types);
void res_function(double t, double *dynamicVars, double *yprime, double *errors, double *constants);
void jac_function(double t, double *dynamicVars, double *yprime, double *pd, double cj, double *constants);
void ic_function(double *dynamicVars, double *constants);
                  
double root(double n,double x){
return pow(x, 1.0/n);
}

double root_0(double n,double x){
return -(log(x)*pow(x, 1.0/n)*1.0/pow(n, 2.0));
}

double root_1(double n,double x){
return pow(x, 1.0/n - 1.0)/n;
}

double cot(double x){
return 1.0/tan(x);
}

double cot_0(double x){
return -(1.0/(pow(cos(x), 2.0)*pow(tan(x), 2.0)));
}

double arccot(double x){
return atan(1.0/x);
}

double arccot_0(double x){
return -(1.0/pow(x, 2.0)/(pow(1.0/x, 2.0) + 1.0));
}

double coth(double x){
return 1.0/tanh(x);
}

double coth_0(double x){
return -(1.0/(pow(cosh(x), 2.0)*pow(tanh(x), 2.0)));
}

double csc(double x){
return 1.0/sin(x);
}

double csc_0(double x){
return -(cos(x)/pow(sin(x), 2.0));
}

double arccsc(double x){
return asin(1.0/x);
}

double arccsc_0(double x){
return -(1.0/pow(x, 2.0)/sqrt(1.0 - pow(1.0/x, 2.0)));
}

double csch(double x){
return 1.0/sinh(x);
}

double csch_0(double x){
return -(cosh(x)/pow(sinh(x), 2.0));
}

double sec(double x){
return 1.0/cos(x);
}

double sec_0(double x){
return sin(x)/pow(cos(x), 2.0);
}

double arcsec(double x){
return acos(1.0/x);
}

double arcsec_0(double x){
return -(1.0/pow(x, 2.0)/sqrt(1.0 - pow(1.0/x, 2.0)));
}

double sech(double x){
return 1.0/cosh(x);
}

double sech_0(double x){
return -(sinh(x)/pow(cosh(x), 2.0));
}


void var_types(int *types) {
types[0] = 1;
types[1] = 1;
types[2] = 1;
types[3] = 1;
types[4] = 1;
types[5] = 1;
types[6] = 1;
types[7] = 1;
types[8] = 1;
types[9] = 1;
types[10] = 1;
types[11] = 1;
types[12] = 1;
types[13] = 1;
types[14] = 1;
}

void res_function(double t, double *dynamicVars, double *yprime, double *errors, double *constants) {
double cell = constants[0];
double RasGapActive = constants[1];
double RapGapActive = constants[2];
double PP2AActive = constants[3];
double Raf1PPtase = constants[4];
double EGF_IC = constants[5];
double NGF_IC = constants[6];
double Total_EGFReceptor = constants[7];
double Total_NGFReceptor = constants[8];
double Total_Sos = constants[9];
double Total_P90Rsk = constants[10];
double Total_Ras = constants[11];
double Total_Raf1 = constants[12];
double Total_BRaf = constants[13];
double Total_Mek = constants[14];
double Total_Erk = constants[15];
double Total_PI3K = constants[16];
double Total_Akt = constants[17];
double Total_C3G = constants[18];
double Total_Rap1 = constants[19];
double Zero = constants[20];
double krbEGF = constants[21];
double kruEGF = constants[22];
double krbNGF = constants[23];
double kruNGF = constants[24];
double kEGF = constants[25];
double KmEGF = constants[26];
double kNGF = constants[27];
double KmNGF = constants[28];
double kdSos = constants[29];
double KmdSos = constants[30];
double kSos = constants[31];
double KmSos = constants[32];
double kRasGap = constants[33];
double KmRasGap = constants[34];
double kRasToRaf1 = constants[35];
double KmRasToRaf1 = constants[36];
double kpRaf1 = constants[37];
double KmpRaf1 = constants[38];
double kpBRaf = constants[39];
double KmpBRaf = constants[40];
double kdMek = constants[41];
double KmdMek = constants[42];
double kpMekCytoplasmic = constants[43];
double KmpMekCytoplasmic = constants[44];
double kdErk = constants[45];
double KmdErk = constants[46];
double kpP90Rsk = constants[47];
double KmpP90Rsk = constants[48];
double kPI3K = constants[49];
double KmPI3K = constants[50];
double kPI3KRas = constants[51];
double KmPI3KRas = constants[52];
double kAkt = constants[53];
double KmAkt = constants[54];
double kdRaf1ByAkt = constants[55];
double KmRaf1ByAkt = constants[56];
double kC3GNGF = constants[57];
double KmC3GNGF = constants[58];
double kC3G = constants[59];
double KmC3G = constants[60];
double kRapGap = constants[61];
double KmRapGap = constants[62];
double kRap1ToBRaf = constants[63];
double KmRap1ToBRaf = constants[64];
double kdRaf1 = constants[65];
double KmdRaf1 = constants[66];
double kdBRaf = constants[67];
double KmdBRaf = constants[68];

double EGF = dynamicVars[0];
double NGF = dynamicVars[1];
double boundEGFReceptor = dynamicVars[2];
double boundNGFReceptor = dynamicVars[3];
double SosActive = dynamicVars[4];
double P90RskActive = dynamicVars[5];
double RasActive = dynamicVars[6];
double Raf1Active = dynamicVars[7];
double BRafActive = dynamicVars[8];
double MekActive = dynamicVars[9];
double ErkActive = dynamicVars[10];
double PI3KActive = dynamicVars[11];
double AktActive = dynamicVars[12];
double C3GActive = dynamicVars[13];
double Rap1Active = dynamicVars[14];
double EGF_prime = yprime[0];
double NGF_prime = yprime[1];
double boundEGFReceptor_prime = yprime[2];
double boundNGFReceptor_prime = yprime[3];
double SosActive_prime = yprime[4];
double P90RskActive_prime = yprime[5];
double RasActive_prime = yprime[6];
double Raf1Active_prime = yprime[7];
double BRafActive_prime = yprime[8];
double MekActive_prime = yprime[9];
double ErkActive_prime = yprime[10];
double PI3KActive_prime = yprime[11];
double AktActive_prime = yprime[12];
double C3GActive_prime = yprime[13];
double Rap1Active_prime = yprime[14];

errors[0] = -EGF*cell*krbEGF*(Total_EGFReceptor - boundEGFReceptor) - EGF_prime + boundEGFReceptor*cell*kruEGF;
errors[1] = -NGF*cell*krbNGF*(Total_NGFReceptor - boundNGFReceptor) - NGF_prime + boundNGFReceptor*cell*kruNGF;
errors[2] = EGF*cell*krbEGF*(Total_EGFReceptor - boundEGFReceptor) - boundEGFReceptor*cell*kruEGF - boundEGFReceptor_prime;
errors[3] = NGF*cell*krbNGF*(Total_NGFReceptor - boundNGFReceptor) - boundNGFReceptor*cell*kruNGF - boundNGFReceptor_prime;
errors[4] = -P90RskActive*SosActive*cell*kdSos/(KmdSos + SosActive) - SosActive_prime + boundEGFReceptor*cell*kEGF*(-SosActive + Total_Sos)/(KmEGF - SosActive + Total_Sos) + boundNGFReceptor*cell*kNGF*(-SosActive + Total_Sos)/(KmNGF - SosActive + Total_Sos);
errors[5] = ErkActive*cell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk) - P90RskActive_prime;
errors[6] = -RasActive*RasGapActive*cell*kRasGap/(KmRasGap + RasActive) - RasActive_prime + SosActive*cell*kSos*(-RasActive + Total_Ras)/(KmSos - RasActive + Total_Ras);
errors[7] = -AktActive*Raf1Active*cell*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) - Raf1Active*Raf1PPtase*cell*kdRaf1/(KmdRaf1 + Raf1Active) - Raf1Active_prime + RasActive*cell*kRasToRaf1*(-Raf1Active + Total_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1);
errors[8] = -BRafActive*Raf1PPtase*cell*kdBRaf/(BRafActive + KmdBRaf) - BRafActive_prime + Rap1Active*cell*kRap1ToBRaf*(-BRafActive + Total_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf);
errors[9] = BRafActive*cell*kpBRaf*(-MekActive + Total_Mek)/(KmpBRaf - MekActive + Total_Mek) - MekActive*PP2AActive*cell*kdMek/(KmdMek + MekActive) - MekActive_prime + Raf1Active*cell*kpRaf1*(-MekActive + Total_Mek)/(KmpRaf1 - MekActive + Total_Mek);
errors[10] = -ErkActive*PP2AActive*cell*kdErk/(ErkActive + KmdErk) - ErkActive_prime + MekActive*cell*kpMekCytoplasmic*(-ErkActive + Total_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk);
errors[11] = -PI3KActive_prime + RasActive*cell*kPI3KRas*(-PI3KActive + Total_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K) + boundEGFReceptor*cell*kPI3K*(-PI3KActive + Total_PI3K)/(KmPI3K - PI3KActive + Total_PI3K);
errors[12] = -AktActive_prime + PI3KActive*cell*kAkt*(-AktActive + Total_Akt)/(-AktActive + KmAkt + Total_Akt);
errors[13] = -C3GActive_prime + boundNGFReceptor*cell*kC3GNGF*(-C3GActive + Total_C3G)/(-C3GActive + KmC3GNGF + Total_C3G);
errors[14] = C3GActive*cell*kC3G*(-Rap1Active + Total_Rap1)/(KmC3G - Rap1Active + Total_Rap1) - Rap1Active*RapGapActive*cell*kRapGap/(KmRapGap + Rap1Active) - Rap1Active_prime;
}

void jac_function(double t, double *dynamicVars, double *yprime, double *pd, double cj, double *constants) {
double cell = constants[0];
double RasGapActive = constants[1];
double RapGapActive = constants[2];
double PP2AActive = constants[3];
double Raf1PPtase = constants[4];
double EGF_IC = constants[5];
double NGF_IC = constants[6];
double Total_EGFReceptor = constants[7];
double Total_NGFReceptor = constants[8];
double Total_Sos = constants[9];
double Total_P90Rsk = constants[10];
double Total_Ras = constants[11];
double Total_Raf1 = constants[12];
double Total_BRaf = constants[13];
double Total_Mek = constants[14];
double Total_Erk = constants[15];
double Total_PI3K = constants[16];
double Total_Akt = constants[17];
double Total_C3G = constants[18];
double Total_Rap1 = constants[19];
double Zero = constants[20];
double krbEGF = constants[21];
double kruEGF = constants[22];
double krbNGF = constants[23];
double kruNGF = constants[24];
double kEGF = constants[25];
double KmEGF = constants[26];
double kNGF = constants[27];
double KmNGF = constants[28];
double kdSos = constants[29];
double KmdSos = constants[30];
double kSos = constants[31];
double KmSos = constants[32];
double kRasGap = constants[33];
double KmRasGap = constants[34];
double kRasToRaf1 = constants[35];
double KmRasToRaf1 = constants[36];
double kpRaf1 = constants[37];
double KmpRaf1 = constants[38];
double kpBRaf = constants[39];
double KmpBRaf = constants[40];
double kdMek = constants[41];
double KmdMek = constants[42];
double kpMekCytoplasmic = constants[43];
double KmpMekCytoplasmic = constants[44];
double kdErk = constants[45];
double KmdErk = constants[46];
double kpP90Rsk = constants[47];
double KmpP90Rsk = constants[48];
double kPI3K = constants[49];
double KmPI3K = constants[50];
double kPI3KRas = constants[51];
double KmPI3KRas = constants[52];
double kAkt = constants[53];
double KmAkt = constants[54];
double kdRaf1ByAkt = constants[55];
double KmRaf1ByAkt = constants[56];
double kC3GNGF = constants[57];
double KmC3GNGF = constants[58];
double kC3G = constants[59];
double KmC3G = constants[60];
double kRapGap = constants[61];
double KmRapGap = constants[62];
double kRap1ToBRaf = constants[63];
double KmRap1ToBRaf = constants[64];
double kdRaf1 = constants[65];
double KmdRaf1 = constants[66];
double kdBRaf = constants[67];
double KmdBRaf = constants[68];

double EGF = dynamicVars[0];
double NGF = dynamicVars[1];
double boundEGFReceptor = dynamicVars[2];
double boundNGFReceptor = dynamicVars[3];
double SosActive = dynamicVars[4];
double P90RskActive = dynamicVars[5];
double RasActive = dynamicVars[6];
double Raf1Active = dynamicVars[7];
double BRafActive = dynamicVars[8];
double MekActive = dynamicVars[9];
double ErkActive = dynamicVars[10];
double PI3KActive = dynamicVars[11];
double AktActive = dynamicVars[12];
double C3GActive = dynamicVars[13];
double Rap1Active = dynamicVars[14];

double EGF_prime = yprime[0];
double NGF_prime = yprime[1];
double boundEGFReceptor_prime = yprime[2];
double boundNGFReceptor_prime = yprime[3];
double SosActive_prime = yprime[4];
double P90RskActive_prime = yprime[5];
double RasActive_prime = yprime[6];
double Raf1Active_prime = yprime[7];
double BRafActive_prime = yprime[8];
double MekActive_prime = yprime[9];
double ErkActive_prime = yprime[10];
double PI3KActive_prime = yprime[11];
double AktActive_prime = yprime[12];
double C3GActive_prime = yprime[13];
double Rap1Active_prime = yprime[14];

memset(pd, 0, sizeof(pd)*225);
pd[0] = (-cell*krbEGF*(Total_EGFReceptor - boundEGFReceptor)) + cj*(-1);
pd[2] = EGF*cell*krbEGF + cell*kruEGF;
pd[16] = (-cell*krbNGF*(Total_NGFReceptor - boundNGFReceptor)) + cj*(-1);
pd[18] = NGF*cell*krbNGF + cell*kruNGF;
pd[30] = cell*krbEGF*(Total_EGFReceptor - boundEGFReceptor);
pd[32] = (-EGF*cell*krbEGF - cell*kruEGF) + cj*(-1);
pd[46] = cell*krbNGF*(Total_NGFReceptor - boundNGFReceptor);
pd[48] = (-NGF*cell*krbNGF - cell*kruNGF) + cj*(-1);
pd[62] = cell*kEGF*(-SosActive + Total_Sos)/(KmEGF - SosActive + Total_Sos);
pd[63] = cell*kNGF*(-SosActive + Total_Sos)/(KmNGF - SosActive + Total_Sos);
pd[64] = (P90RskActive*SosActive*cell*kdSos/pow(KmdSos + SosActive, 2) - P90RskActive*cell*kdSos/(KmdSos + SosActive) + boundEGFReceptor*cell*kEGF*(-SosActive + Total_Sos)/pow(KmEGF - SosActive + Total_Sos, 2) - boundEGFReceptor*cell*kEGF/(KmEGF - SosActive + Total_Sos) + boundNGFReceptor*cell*kNGF*(-SosActive + Total_Sos)/pow(KmNGF - SosActive + Total_Sos, 2) - boundNGFReceptor*cell*kNGF/(KmNGF - SosActive + Total_Sos)) + cj*(-1);
pd[65] = -SosActive*cell*kdSos/(KmdSos + SosActive);
pd[80] = (ErkActive*cell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 2) - ErkActive*cell*kpP90Rsk/(KmpP90Rsk - P90RskActive + Total_P90Rsk)) + cj*(-1);
pd[85] = cell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk);
pd[94] = cell*kSos*(-RasActive + Total_Ras)/(KmSos - RasActive + Total_Ras);
pd[96] = (RasActive*RasGapActive*cell*kRasGap/pow(KmRasGap + RasActive, 2) - RasGapActive*cell*kRasGap/(KmRasGap + RasActive) + SosActive*cell*kSos*(-RasActive + Total_Ras)/pow(KmSos - RasActive + Total_Ras, 2) - SosActive*cell*kSos/(KmSos - RasActive + Total_Ras)) + cj*(-1);
pd[111] = cell*kRasToRaf1*(-Raf1Active + Total_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1);
pd[112] = (AktActive*Raf1Active*cell*kdRaf1ByAkt/pow(KmRaf1ByAkt + Raf1Active, 2) - AktActive*cell*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) + Raf1Active*Raf1PPtase*cell*kdRaf1/pow(KmdRaf1 + Raf1Active, 2) - Raf1PPtase*cell*kdRaf1/(KmdRaf1 + Raf1Active) + RasActive*cell*kRasToRaf1*(-Raf1Active + Total_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 2) - RasActive*cell*kRasToRaf1/(KmRasToRaf1 - Raf1Active + Total_Raf1)) + cj*(-1);
pd[117] = -Raf1Active*cell*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active);
pd[128] = (BRafActive*Raf1PPtase*cell*kdBRaf/pow(BRafActive + KmdBRaf, 2) - Raf1PPtase*cell*kdBRaf/(BRafActive + KmdBRaf) + Rap1Active*cell*kRap1ToBRaf*(-BRafActive + Total_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 2) - Rap1Active*cell*kRap1ToBRaf/(-BRafActive + KmRap1ToBRaf + Total_BRaf)) + cj*(-1);
pd[134] = cell*kRap1ToBRaf*(-BRafActive + Total_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf);
pd[142] = cell*kpRaf1*(-MekActive + Total_Mek)/(KmpRaf1 - MekActive + Total_Mek);
pd[143] = cell*kpBRaf*(-MekActive + Total_Mek)/(KmpBRaf - MekActive + Total_Mek);
pd[144] = (BRafActive*cell*kpBRaf*(-MekActive + Total_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 2) - BRafActive*cell*kpBRaf/(KmpBRaf - MekActive + Total_Mek) + MekActive*PP2AActive*cell*kdMek/pow(KmdMek + MekActive, 2) - PP2AActive*cell*kdMek/(KmdMek + MekActive) + Raf1Active*cell*kpRaf1*(-MekActive + Total_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 2) - Raf1Active*cell*kpRaf1/(KmpRaf1 - MekActive + Total_Mek)) + cj*(-1);
pd[159] = cell*kpMekCytoplasmic*(-ErkActive + Total_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk);
pd[160] = (ErkActive*PP2AActive*cell*kdErk/pow(ErkActive + KmdErk, 2) + MekActive*cell*kpMekCytoplasmic*(-ErkActive + Total_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 2) - MekActive*cell*kpMekCytoplasmic/(-ErkActive + KmpMekCytoplasmic + Total_Erk) - PP2AActive*cell*kdErk/(ErkActive + KmdErk)) + cj*(-1);
pd[167] = cell*kPI3K*(-PI3KActive + Total_PI3K)/(KmPI3K - PI3KActive + Total_PI3K);
pd[171] = cell*kPI3KRas*(-PI3KActive + Total_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K);
pd[176] = (RasActive*cell*kPI3KRas*(-PI3KActive + Total_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 2) - RasActive*cell*kPI3KRas/(KmPI3KRas - PI3KActive + Total_PI3K) + boundEGFReceptor*cell*kPI3K*(-PI3KActive + Total_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 2) - boundEGFReceptor*cell*kPI3K/(KmPI3K - PI3KActive + Total_PI3K)) + cj*(-1);
pd[191] = cell*kAkt*(-AktActive + Total_Akt)/(-AktActive + KmAkt + Total_Akt);
pd[192] = (PI3KActive*cell*kAkt*(-AktActive + Total_Akt)/pow(-AktActive + KmAkt + Total_Akt, 2) - PI3KActive*cell*kAkt/(-AktActive + KmAkt + Total_Akt)) + cj*(-1);
pd[198] = cell*kC3GNGF*(-C3GActive + Total_C3G)/(-C3GActive + KmC3GNGF + Total_C3G);
pd[208] = (boundNGFReceptor*cell*kC3GNGF*(-C3GActive + Total_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 2) - boundNGFReceptor*cell*kC3GNGF/(-C3GActive + KmC3GNGF + Total_C3G)) + cj*(-1);
pd[223] = cell*kC3G*(-Rap1Active + Total_Rap1)/(KmC3G - Rap1Active + Total_Rap1);
pd[224] = (C3GActive*cell*kC3G*(-Rap1Active + Total_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 2) - C3GActive*cell*kC3G/(KmC3G - Rap1Active + Total_Rap1) + Rap1Active*RapGapActive*cell*kRapGap/pow(KmRapGap + Rap1Active, 2) - RapGapActive*cell*kRapGap/(KmRapGap + Rap1Active)) + cj*(-1);
}
void ic_function(double *dynamicVars, double *constants) {
double cell = constants[0];
double RasGapActive = constants[1];
double RapGapActive = constants[2];
double PP2AActive = constants[3];
double Raf1PPtase = constants[4];
double EGF_IC = constants[5];
double NGF_IC = constants[6];
double Total_EGFReceptor = constants[7];
double Total_NGFReceptor = constants[8];
double Total_Sos = constants[9];
double Total_P90Rsk = constants[10];
double Total_Ras = constants[11];
double Total_Raf1 = constants[12];
double Total_BRaf = constants[13];
double Total_Mek = constants[14];
double Total_Erk = constants[15];
double Total_PI3K = constants[16];
double Total_Akt = constants[17];
double Total_C3G = constants[18];
double Total_Rap1 = constants[19];
double Zero = constants[20];
double krbEGF = constants[21];
double kruEGF = constants[22];
double krbNGF = constants[23];
double kruNGF = constants[24];
double kEGF = constants[25];
double KmEGF = constants[26];
double kNGF = constants[27];
double KmNGF = constants[28];
double kdSos = constants[29];
double KmdSos = constants[30];
double kSos = constants[31];
double KmSos = constants[32];
double kRasGap = constants[33];
double KmRasGap = constants[34];
double kRasToRaf1 = constants[35];
double KmRasToRaf1 = constants[36];
double kpRaf1 = constants[37];
double KmpRaf1 = constants[38];
double kpBRaf = constants[39];
double KmpBRaf = constants[40];
double kdMek = constants[41];
double KmdMek = constants[42];
double kpMekCytoplasmic = constants[43];
double KmpMekCytoplasmic = constants[44];
double kdErk = constants[45];
double KmdErk = constants[46];
double kpP90Rsk = constants[47];
double KmpP90Rsk = constants[48];
double kPI3K = constants[49];
double KmPI3K = constants[50];
double kPI3KRas = constants[51];
double KmPI3KRas = constants[52];
double kAkt = constants[53];
double KmAkt = constants[54];
double kdRaf1ByAkt = constants[55];
double KmRaf1ByAkt = constants[56];
double kC3GNGF = constants[57];
double KmC3GNGF = constants[58];
double kC3G = constants[59];
double KmC3G = constants[60];
double kRapGap = constants[61];
double KmRapGap = constants[62];
double kRap1ToBRaf = constants[63];
double KmRap1ToBRaf = constants[64];
double kdRaf1 = constants[65];
double KmdRaf1 = constants[66];
double kdBRaf = constants[67];
double KmdBRaf = constants[68];

dynamicVars[0] = EGF_IC;
dynamicVars[1] = NGF_IC;
dynamicVars[2] = Zero;
dynamicVars[3] = Zero;
dynamicVars[4] = Zero;
dynamicVars[5] = Zero;
dynamicVars[6] = Zero;
dynamicVars[7] = Zero;
dynamicVars[8] = Zero;
dynamicVars[9] = Zero;
dynamicVars[10] = Zero;
dynamicVars[11] = Zero;
dynamicVars[12] = Zero;
dynamicVars[13] = Zero;
dynamicVars[14] = Zero;
}
