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
types[15] = 1;
types[16] = 1;
types[17] = 1;
types[18] = 1;
types[19] = 1;
types[20] = 1;
types[21] = 1;
types[22] = 1;
types[23] = 1;
types[24] = 1;
types[25] = 1;
types[26] = 1;
types[27] = 1;
types[28] = 1;
types[29] = 1;
types[30] = 1;
types[31] = 1;
types[32] = 1;
types[33] = 1;
types[34] = 1;
types[35] = 1;
types[36] = 1;
types[37] = 1;
types[38] = 1;
types[39] = 1;
types[40] = 1;
types[41] = 1;
types[42] = 1;
types[43] = 1;
types[44] = 1;
types[45] = 1;
types[46] = 1;
types[47] = 1;
types[48] = 1;
types[49] = 1;
types[50] = 1;
types[51] = 1;
types[52] = 1;
types[53] = 1;
}

void res_function(double t, double *dynamicVars, double *yprime, double *errors, double *constants) {
double cell = constants[0];
double RasGapTotal = constants[1];
double RapGapTotal = constants[2];
double PP2ATotal = constants[3];
double Raf1PPtaseTotal = constants[4];
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
double k00f = constants[21];
double k00b = constants[22];
double k01f = constants[23];
double k01b = constants[24];
double k02f = constants[25];
double k02b = constants[26];
double k03f = constants[27];
double k04f = constants[28];
double k04b = constants[29];
double k05f = constants[30];
double k06f = constants[31];
double k06b = constants[32];
double k07f = constants[33];
double k08f = constants[34];
double k08b = constants[35];
double k09f = constants[36];
double k10f = constants[37];
double k10b = constants[38];
double k11f = constants[39];
double k12f = constants[40];
double k12b = constants[41];
double k13f = constants[42];
double k14f = constants[43];
double k14b = constants[44];
double k15f = constants[45];
double k16f = constants[46];
double k16b = constants[47];
double k17f = constants[48];
double k18f = constants[49];
double k18b = constants[50];
double k19f = constants[51];
double k20f = constants[52];
double k20b = constants[53];
double k21f = constants[54];
double k22f = constants[55];
double k22b = constants[56];
double k23f = constants[57];
double k24f = constants[58];
double k24b = constants[59];
double k25f = constants[60];
double k26f = constants[61];
double k26b = constants[62];
double k27f = constants[63];
double k28f = constants[64];
double k28b = constants[65];
double k29f = constants[66];
double k30f = constants[67];
double k30b = constants[68];
double k31f = constants[69];
double k32f = constants[70];
double k32b = constants[71];
double k33f = constants[72];
double k34f = constants[73];
double k34b = constants[74];
double k35f = constants[75];
double k36f = constants[76];
double k36b = constants[77];
double k37f = constants[78];
double k38f = constants[79];
double k38b = constants[80];
double k39f = constants[81];
double k40f = constants[82];
double k40b = constants[83];
double k41f = constants[84];
double k42f = constants[85];
double k42b = constants[86];
double k43f = constants[87];
double k44f = constants[88];
double k44b = constants[89];
double k45f = constants[90];

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
double freeEGFReceptor = dynamicVars[15];
double freeNGFReceptor = dynamicVars[16];
double SosInactive = dynamicVars[17];
double P90RskInactive = dynamicVars[18];
double RasInactive = dynamicVars[19];
double Raf1Inactive = dynamicVars[20];
double BRafInactive = dynamicVars[21];
double MekInactive = dynamicVars[22];
double ErkInactive = dynamicVars[23];
double PI3KInactive = dynamicVars[24];
double AktInactive = dynamicVars[25];
double C3GInactive = dynamicVars[26];
double Rap1Inactive = dynamicVars[27];
double RasPI3K = dynamicVars[28];
double Raf1PPtaseRaf1 = dynamicVars[29];
double BRafMek = dynamicVars[30];
double boundNGFReceptorC3G = dynamicVars[31];
double P90RskSos = dynamicVars[32];
double Raf1PPtaseBRaf = dynamicVars[33];
double PP2AActive = dynamicVars[34];
double RasGapRas = dynamicVars[35];
double ErkP90Rsk = dynamicVars[36];
double RasRaf1 = dynamicVars[37];
double MekErk = dynamicVars[38];
double RapGapRap1 = dynamicVars[39];
double RapGapActive = dynamicVars[40];
double AktRaf1 = dynamicVars[41];
double Raf1PPtase = dynamicVars[42];
double C3GRap1 = dynamicVars[43];
double Raf1Mek = dynamicVars[44];
double SosRas = dynamicVars[45];
double boundEGFReceptorPI3K = dynamicVars[46];
double RasGapActive = dynamicVars[47];
double boundEGFReceptorSos = dynamicVars[48];
double boundNGFReceptorSos = dynamicVars[49];
double Rap1BRaf = dynamicVars[50];
double PP2AMek = dynamicVars[51];
double PP2AErk = dynamicVars[52];
double PI3KAkt = dynamicVars[53];
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
double freeEGFReceptor_prime = yprime[15];
double freeNGFReceptor_prime = yprime[16];
double SosInactive_prime = yprime[17];
double P90RskInactive_prime = yprime[18];
double RasInactive_prime = yprime[19];
double Raf1Inactive_prime = yprime[20];
double BRafInactive_prime = yprime[21];
double MekInactive_prime = yprime[22];
double ErkInactive_prime = yprime[23];
double PI3KInactive_prime = yprime[24];
double AktInactive_prime = yprime[25];
double C3GInactive_prime = yprime[26];
double Rap1Inactive_prime = yprime[27];
double RasPI3K_prime = yprime[28];
double Raf1PPtaseRaf1_prime = yprime[29];
double BRafMek_prime = yprime[30];
double boundNGFReceptorC3G_prime = yprime[31];
double P90RskSos_prime = yprime[32];
double Raf1PPtaseBRaf_prime = yprime[33];
double PP2AActive_prime = yprime[34];
double RasGapRas_prime = yprime[35];
double ErkP90Rsk_prime = yprime[36];
double RasRaf1_prime = yprime[37];
double MekErk_prime = yprime[38];
double RapGapRap1_prime = yprime[39];
double RapGapActive_prime = yprime[40];
double AktRaf1_prime = yprime[41];
double Raf1PPtase_prime = yprime[42];
double C3GRap1_prime = yprime[43];
double Raf1Mek_prime = yprime[44];
double SosRas_prime = yprime[45];
double boundEGFReceptorPI3K_prime = yprime[46];
double RasGapActive_prime = yprime[47];
double boundEGFReceptorSos_prime = yprime[48];
double boundNGFReceptorSos_prime = yprime[49];
double Rap1BRaf_prime = yprime[50];
double PP2AMek_prime = yprime[51];
double PP2AErk_prime = yprime[52];
double PI3KAkt_prime = yprime[53];

errors[0] = -C3GInactive*boundNGFReceptor*k40f + NGF*freeNGFReceptor*k01f - SosInactive*boundNGFReceptor*k04f - boundNGFReceptor*k01b + boundNGFReceptorC3G*k40b + boundNGFReceptorC3G*k41f + boundNGFReceptorSos*k04b + boundNGFReceptorSos*k05f - boundNGFReceptor_prime;
errors[1] = -NGF*freeNGFReceptor*k01f - NGF_prime + boundNGFReceptor*k01b;
errors[2] = PI3KInactive*RasActive*k36f - RasPI3K*k36b - RasPI3K*k37f - RasPI3K_prime;
errors[3] = -ErkInactive*MekActive*k30f - ErkInactive_prime + MekErk*k30b + PP2AErk*k33f;
errors[4] = BRafMek*k27f - ErkInactive*MekActive*k30f - MekActive*PP2AActive*k28f - MekActive_prime + MekErk*k30b + MekErk*k31f + PP2AMek*k28b + Raf1Mek*k25f;
errors[5] = SosInactive*boundNGFReceptor*k04f - boundNGFReceptorSos*k04b - boundNGFReceptorSos*k05f - boundNGFReceptorSos_prime;
errors[6] = Raf1Active*Raf1PPtase*k16f - Raf1PPtaseRaf1*k16b - Raf1PPtaseRaf1*k17f - Raf1PPtaseRaf1_prime;
errors[7] = BRafActive*MekInactive*k26f - BRafMek*k26b - BRafMek*k27f - BRafMek_prime;
errors[8] = C3GInactive*boundNGFReceptor*k40f - boundNGFReceptorC3G*k40b - boundNGFReceptorC3G*k41f - boundNGFReceptorC3G_prime;
errors[9] = P90RskActive*SosActive*k06f - P90RskSos*k06b - P90RskSos*k07f - P90RskSos_prime;
errors[10] = -ErkActive*P90RskInactive*k08f + ErkP90Rsk*k08b - P90RskInactive_prime;
errors[11] = -EGF*freeEGFReceptor*k00f - EGF_prime + boundEGFReceptor*k00b;
errors[12] = BRafActive*Raf1PPtase*k22f - Raf1PPtaseBRaf*k22b - Raf1PPtaseBRaf*k23f - Raf1PPtaseBRaf_prime;
errors[13] = -ErkActive*PP2AActive*k32f - MekActive*PP2AActive*k28f - PP2AActive_prime + PP2AErk*k32b + PP2AErk*k33f + PP2AMek*k28b + PP2AMek*k29f;
errors[14] = -AktInactive*PI3KActive*k38f - AktInactive_prime + PI3KAkt*k38b;
errors[15] = -AktActive*Raf1Active*k18f - AktActive_prime + AktRaf1*k18b + AktRaf1*k19f + PI3KAkt*k39f;
errors[16] = RasActive*RasGapActive*k12f - RasGapRas*k12b - RasGapRas*k13f - RasGapRas_prime;
errors[17] = -P90RskActive*SosActive*k06f + P90RskSos*k06b - RasInactive*SosActive*k10f - SosActive_prime + SosRas*k10b + SosRas*k11f + boundEGFReceptorSos*k03f + boundNGFReceptorSos*k05f;
errors[18] = ErkActive*P90RskInactive*k08f - ErkP90Rsk*k08b - ErkP90Rsk*k09f - ErkP90Rsk_prime;
errors[19] = Raf1Inactive*RasActive*k14f - RasRaf1*k14b - RasRaf1*k15f - RasRaf1_prime;
errors[20] = EGF*freeEGFReceptor*k00f - PI3KInactive*boundEGFReceptor*k34f - SosInactive*boundEGFReceptor*k02f - boundEGFReceptor*k00b + boundEGFReceptorPI3K*k34b + boundEGFReceptorPI3K*k35f + boundEGFReceptorSos*k02b + boundEGFReceptorSos*k03f - boundEGFReceptor_prime;
errors[21] = -AktActive*Raf1Active*k18f + AktRaf1*k18b - MekInactive*Raf1Active*k24f - Raf1Active*Raf1PPtase*k16f - Raf1Active_prime + Raf1Mek*k24b + Raf1Mek*k25f + Raf1PPtaseRaf1*k16b + RasRaf1*k15f;
errors[22] = ErkInactive*MekActive*k30f - MekErk*k30b - MekErk*k31f - MekErk_prime;
errors[23] = Rap1Active*RapGapActive*k44f - RapGapRap1*k44b - RapGapRap1*k45f - RapGapRap1_prime;
errors[24] = P90RskSos*k07f - SosInactive*boundEGFReceptor*k02f - SosInactive*boundNGFReceptor*k04f - SosInactive_prime + boundEGFReceptorSos*k02b + boundNGFReceptorSos*k04b;
errors[25] = -Rap1Active*RapGapActive*k44f - RapGapActive_prime + RapGapRap1*k44b + RapGapRap1*k45f;
errors[26] = AktActive*Raf1Active*k18f - AktRaf1*k18b - AktRaf1*k19f - AktRaf1_prime;
errors[27] = -EGF*freeEGFReceptor*k00f + boundEGFReceptor*k00b - freeEGFReceptor_prime;
errors[28] = -BRafActive*Raf1PPtase*k22f - Raf1Active*Raf1PPtase*k16f + Raf1PPtaseBRaf*k22b + Raf1PPtaseBRaf*k23f + Raf1PPtaseRaf1*k16b + Raf1PPtaseRaf1*k17f - Raf1PPtase_prime;
errors[29] = C3GActive*Rap1Inactive*k42f - C3GRap1*k42b - C3GRap1*k43f - C3GRap1_prime;
errors[30] = -BRafInactive*Rap1Active*k20f + C3GRap1*k43f - Rap1Active*RapGapActive*k44f - Rap1Active_prime + Rap1BRaf*k20b + Rap1BRaf*k21f + RapGapRap1*k44b;
errors[31] = -BRafActive*MekInactive*k26f + BRafMek*k26b - MekInactive*Raf1Active*k24f - MekInactive_prime + PP2AMek*k29f + Raf1Mek*k24b;
errors[32] = MekInactive*Raf1Active*k24f - Raf1Mek*k24b - Raf1Mek*k25f - Raf1Mek_prime;
errors[33] = RasInactive*SosActive*k10f - SosRas*k10b - SosRas*k11f - SosRas_prime;
errors[34] = ErkP90Rsk*k09f - P90RskActive*SosActive*k06f - P90RskActive_prime + P90RskSos*k06b + P90RskSos*k07f;
errors[35] = ErkActive*PP2AActive*k32f - PP2AErk*k32b - PP2AErk*k33f - PP2AErk_prime;
errors[36] = PI3KInactive*boundEGFReceptor*k34f - boundEGFReceptorPI3K*k34b - boundEGFReceptorPI3K*k35f - boundEGFReceptorPI3K_prime;
errors[37] = -AktInactive*PI3KActive*k38f - PI3KActive_prime + PI3KAkt*k38b + PI3KAkt*k39f + RasPI3K*k37f + boundEGFReceptorPI3K*k35f;
errors[38] = -RasActive*RasGapActive*k12f - RasGapActive_prime + RasGapRas*k12b + RasGapRas*k13f;
errors[39] = -BRafActive*MekInactive*k26f - BRafActive*Raf1PPtase*k22f - BRafActive_prime + BRafMek*k26b + BRafMek*k27f + Raf1PPtaseBRaf*k22b + Rap1BRaf*k21f;
errors[40] = SosInactive*boundEGFReceptor*k02f - boundEGFReceptorSos*k02b - boundEGFReceptorSos*k03f - boundEGFReceptorSos_prime;
errors[41] = -C3GActive*Rap1Inactive*k42f + C3GRap1*k42b - Rap1Inactive_prime + RapGapRap1*k45f;
errors[42] = -ErkActive*P90RskInactive*k08f - ErkActive*PP2AActive*k32f - ErkActive_prime + ErkP90Rsk*k08b + ErkP90Rsk*k09f + MekErk*k31f + PP2AErk*k32b;
errors[43] = -NGF*freeNGFReceptor*k01f + boundNGFReceptor*k01b - freeNGFReceptor_prime;
errors[44] = -C3GInactive*boundNGFReceptor*k40f - C3GInactive_prime + boundNGFReceptorC3G*k40b;
errors[45] = AktRaf1*k19f - Raf1Inactive*RasActive*k14f - Raf1Inactive_prime + Raf1PPtaseRaf1*k17f + RasRaf1*k14b;
errors[46] = -PI3KInactive*RasActive*k36f - PI3KInactive*boundEGFReceptor*k34f - PI3KInactive_prime + RasPI3K*k36b + boundEGFReceptorPI3K*k34b;
errors[47] = -BRafInactive*Rap1Active*k20f - BRafInactive_prime + Raf1PPtaseBRaf*k23f + Rap1BRaf*k20b;
errors[48] = BRafInactive*Rap1Active*k20f - Rap1BRaf*k20b - Rap1BRaf*k21f - Rap1BRaf_prime;
errors[49] = -PI3KInactive*RasActive*k36f - Raf1Inactive*RasActive*k14f - RasActive*RasGapActive*k12f - RasActive_prime + RasGapRas*k12b + RasPI3K*k36b + RasPI3K*k37f + RasRaf1*k14b + RasRaf1*k15f + SosRas*k11f;
errors[50] = RasGapRas*k13f - RasInactive*SosActive*k10f - RasInactive_prime + SosRas*k10b;
errors[51] = MekActive*PP2AActive*k28f - PP2AMek*k28b - PP2AMek*k29f - PP2AMek_prime;
errors[52] = AktInactive*PI3KActive*k38f - PI3KAkt*k38b - PI3KAkt*k39f - PI3KAkt_prime;
errors[53] = -C3GActive*Rap1Inactive*k42f - C3GActive_prime + C3GRap1*k42b + C3GRap1*k43f + boundNGFReceptorC3G*k41f;
}

void jac_function(double t, double *dynamicVars, double *yprime, double *pd, double cj, double *constants) {
double cell = constants[0];
double RasGapTotal = constants[1];
double RapGapTotal = constants[2];
double PP2ATotal = constants[3];
double Raf1PPtaseTotal = constants[4];
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
double k00f = constants[21];
double k00b = constants[22];
double k01f = constants[23];
double k01b = constants[24];
double k02f = constants[25];
double k02b = constants[26];
double k03f = constants[27];
double k04f = constants[28];
double k04b = constants[29];
double k05f = constants[30];
double k06f = constants[31];
double k06b = constants[32];
double k07f = constants[33];
double k08f = constants[34];
double k08b = constants[35];
double k09f = constants[36];
double k10f = constants[37];
double k10b = constants[38];
double k11f = constants[39];
double k12f = constants[40];
double k12b = constants[41];
double k13f = constants[42];
double k14f = constants[43];
double k14b = constants[44];
double k15f = constants[45];
double k16f = constants[46];
double k16b = constants[47];
double k17f = constants[48];
double k18f = constants[49];
double k18b = constants[50];
double k19f = constants[51];
double k20f = constants[52];
double k20b = constants[53];
double k21f = constants[54];
double k22f = constants[55];
double k22b = constants[56];
double k23f = constants[57];
double k24f = constants[58];
double k24b = constants[59];
double k25f = constants[60];
double k26f = constants[61];
double k26b = constants[62];
double k27f = constants[63];
double k28f = constants[64];
double k28b = constants[65];
double k29f = constants[66];
double k30f = constants[67];
double k30b = constants[68];
double k31f = constants[69];
double k32f = constants[70];
double k32b = constants[71];
double k33f = constants[72];
double k34f = constants[73];
double k34b = constants[74];
double k35f = constants[75];
double k36f = constants[76];
double k36b = constants[77];
double k37f = constants[78];
double k38f = constants[79];
double k38b = constants[80];
double k39f = constants[81];
double k40f = constants[82];
double k40b = constants[83];
double k41f = constants[84];
double k42f = constants[85];
double k42b = constants[86];
double k43f = constants[87];
double k44f = constants[88];
double k44b = constants[89];
double k45f = constants[90];

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
double freeEGFReceptor = dynamicVars[15];
double freeNGFReceptor = dynamicVars[16];
double SosInactive = dynamicVars[17];
double P90RskInactive = dynamicVars[18];
double RasInactive = dynamicVars[19];
double Raf1Inactive = dynamicVars[20];
double BRafInactive = dynamicVars[21];
double MekInactive = dynamicVars[22];
double ErkInactive = dynamicVars[23];
double PI3KInactive = dynamicVars[24];
double AktInactive = dynamicVars[25];
double C3GInactive = dynamicVars[26];
double Rap1Inactive = dynamicVars[27];
double RasPI3K = dynamicVars[28];
double Raf1PPtaseRaf1 = dynamicVars[29];
double BRafMek = dynamicVars[30];
double boundNGFReceptorC3G = dynamicVars[31];
double P90RskSos = dynamicVars[32];
double Raf1PPtaseBRaf = dynamicVars[33];
double PP2AActive = dynamicVars[34];
double RasGapRas = dynamicVars[35];
double ErkP90Rsk = dynamicVars[36];
double RasRaf1 = dynamicVars[37];
double MekErk = dynamicVars[38];
double RapGapRap1 = dynamicVars[39];
double RapGapActive = dynamicVars[40];
double AktRaf1 = dynamicVars[41];
double Raf1PPtase = dynamicVars[42];
double C3GRap1 = dynamicVars[43];
double Raf1Mek = dynamicVars[44];
double SosRas = dynamicVars[45];
double boundEGFReceptorPI3K = dynamicVars[46];
double RasGapActive = dynamicVars[47];
double boundEGFReceptorSos = dynamicVars[48];
double boundNGFReceptorSos = dynamicVars[49];
double Rap1BRaf = dynamicVars[50];
double PP2AMek = dynamicVars[51];
double PP2AErk = dynamicVars[52];
double PI3KAkt = dynamicVars[53];

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
double freeEGFReceptor_prime = yprime[15];
double freeNGFReceptor_prime = yprime[16];
double SosInactive_prime = yprime[17];
double P90RskInactive_prime = yprime[18];
double RasInactive_prime = yprime[19];
double Raf1Inactive_prime = yprime[20];
double BRafInactive_prime = yprime[21];
double MekInactive_prime = yprime[22];
double ErkInactive_prime = yprime[23];
double PI3KInactive_prime = yprime[24];
double AktInactive_prime = yprime[25];
double C3GInactive_prime = yprime[26];
double Rap1Inactive_prime = yprime[27];
double RasPI3K_prime = yprime[28];
double Raf1PPtaseRaf1_prime = yprime[29];
double BRafMek_prime = yprime[30];
double boundNGFReceptorC3G_prime = yprime[31];
double P90RskSos_prime = yprime[32];
double Raf1PPtaseBRaf_prime = yprime[33];
double PP2AActive_prime = yprime[34];
double RasGapRas_prime = yprime[35];
double ErkP90Rsk_prime = yprime[36];
double RasRaf1_prime = yprime[37];
double MekErk_prime = yprime[38];
double RapGapRap1_prime = yprime[39];
double RapGapActive_prime = yprime[40];
double AktRaf1_prime = yprime[41];
double Raf1PPtase_prime = yprime[42];
double C3GRap1_prime = yprime[43];
double Raf1Mek_prime = yprime[44];
double SosRas_prime = yprime[45];
double boundEGFReceptorPI3K_prime = yprime[46];
double RasGapActive_prime = yprime[47];
double boundEGFReceptorSos_prime = yprime[48];
double boundNGFReceptorSos_prime = yprime[49];
double Rap1BRaf_prime = yprime[50];
double PP2AMek_prime = yprime[51];
double PP2AErk_prime = yprime[52];
double PI3KAkt_prime = yprime[53];

memset(pd, 0, sizeof(pd)*2916);
pd[1] = freeNGFReceptor*k01f;
pd[3] = (-C3GInactive*k40f - SosInactive*k04f - k01b) + cj*(-1);
pd[16] = NGF*k01f;
pd[17] = -boundNGFReceptor*k04f;
pd[26] = -boundNGFReceptor*k40f;
pd[31] = k40b + k41f;
pd[49] = k04b + k05f;
pd[55] = (-freeNGFReceptor*k01f) + cj*(-1);
pd[57] = k01b;
pd[70] = -NGF*k01f;
pd[114] = PI3KInactive*k36f;
pd[132] = RasActive*k36f;
pd[136] = (-k36b - k37f) + cj*(-1);
pd[171] = -ErkInactive*k30f;
pd[185] = (-MekActive*k30f) + cj*(-1);
pd[200] = k30b;
pd[214] = k33f;
pd[225] = (-ErkInactive*k30f - PP2AActive*k28f) + cj*(-1);
pd[239] = -MekActive*k30f;
pd[246] = k27f;
pd[250] = -MekActive*k28f;
pd[254] = k30b + k31f;
pd[260] = k25f;
pd[267] = k28b;
pd[273] = SosInactive*k04f;
pd[287] = boundNGFReceptor*k04f;
pd[319] = (-k04b - k05f) + cj*(-1);
pd[331] = Raf1PPtase*k16f;
pd[353] = (-k16b - k17f) + cj*(-1);
pd[366] = Raf1Active*k16f;
pd[386] = MekInactive*k26f;
pd[400] = BRafActive*k26f;
pd[408] = (-k26b - k27f) + cj*(-1);
pd[435] = C3GInactive*k40f;
pd[458] = boundNGFReceptor*k40f;
pd[463] = (-k40b - k41f) + cj*(-1);
pd[490] = P90RskActive*k06f;
pd[491] = SosActive*k06f;
pd[518] = (-k06b - k07f) + cj*(-1);
pd[550] = -P90RskInactive*k08f;
pd[558] = (-ErkActive*k08f) + cj*(-1);
pd[576] = k08b;
pd[594] = (-freeEGFReceptor*k00f) + cj*(-1);
pd[596] = k00b;
pd[609] = -EGF*k00f;
pd[656] = Raf1PPtase*k22f;
pd[681] = (-k22b - k23f) + cj*(-1);
pd[690] = BRafActive*k22f;
pd[711] = -PP2AActive*k28f;
pd[712] = -PP2AActive*k32f;
pd[736] = (-ErkActive*k32f - MekActive*k28f) + cj*(-1);
pd[753] = k28b + k29f;
pd[754] = k32b + k33f;
pd[767] = -AktInactive*k38f;
pd[781] = (-PI3KActive*k38f) + cj*(-1);
pd[809] = k38b;
pd[817] = -AktActive*k18f;
pd[822] = (-Raf1Active*k18f) + cj*(-1);
pd[851] = k18b + k19f;
pd[863] = k39f;
pd[870] = RasGapActive*k12f;
pd[899] = (-k12b - k13f) + cj*(-1);
pd[911] = RasActive*k12f;
pd[922] = (-P90RskActive*k06f - RasInactive*k10f) + cj*(-1);
pd[923] = -SosActive*k06f;
pd[937] = -SosActive*k10f;
pd[950] = k06b;
pd[963] = k10b + k11f;
pd[966] = k03f;
pd[967] = k05f;
pd[982] = P90RskInactive*k08f;
pd[990] = ErkActive*k08f;
pd[1008] = (-k08b - k09f) + cj*(-1);
pd[1032] = Raf1Inactive*k14f;
pd[1046] = RasActive*k14f;
pd[1063] = (-k14b - k15f) + cj*(-1);
pd[1080] = freeEGFReceptor*k00f;
pd[1082] = (-PI3KInactive*k34f - SosInactive*k02f - k00b) + cj*(-1);
pd[1095] = EGF*k00f;
pd[1097] = -boundEGFReceptor*k02f;
pd[1104] = -boundEGFReceptor*k34f;
pd[1126] = k34b + k35f;
pd[1128] = k02b + k03f;
pd[1141] = (-AktActive*k18f - MekInactive*k24f - Raf1PPtase*k16f) + cj*(-1);
pd[1146] = -Raf1Active*k18f;
pd[1156] = -Raf1Active*k24f;
pd[1163] = k16b;
pd[1171] = k15f;
pd[1175] = k18b;
pd[1176] = -Raf1Active*k16f;
pd[1178] = k24b + k25f;
pd[1197] = ErkInactive*k30f;
pd[1211] = MekActive*k30f;
pd[1226] = (-k30b - k31f) + cj*(-1);
pd[1256] = RapGapActive*k44f;
pd[1281] = (-k44b - k45f) + cj*(-1);
pd[1282] = Rap1Active*k44f;
pd[1298] = -SosInactive*k02f;
pd[1299] = -SosInactive*k04f;
pd[1313] = (-boundEGFReceptor*k02f - boundNGFReceptor*k04f) + cj*(-1);
pd[1328] = k07f;
pd[1344] = k02b;
pd[1345] = k04b;
pd[1364] = -RapGapActive*k44f;
pd[1389] = k44b + k45f;
pd[1390] = (-Rap1Active*k44f) + cj*(-1);
pd[1411] = AktActive*k18f;
pd[1416] = Raf1Active*k18f;
pd[1445] = (-k18b - k19f) + cj*(-1);
pd[1458] = -freeEGFReceptor*k00f;
pd[1460] = k00b;
pd[1473] = (-EGF*k00f) + cj*(-1);
pd[1519] = -Raf1PPtase*k16f;
pd[1520] = -Raf1PPtase*k22f;
pd[1541] = k16b + k17f;
pd[1545] = k22b + k23f;
pd[1554] = (-BRafActive*k22f - Raf1Active*k16f) + cj*(-1);
pd[1579] = Rap1Inactive*k42f;
pd[1593] = C3GActive*k42f;
pd[1609] = (-k42b - k43f) + cj*(-1);
pd[1634] = (-BRafInactive*k20f - RapGapActive*k44f) + cj*(-1);
pd[1641] = -Rap1Active*k20f;
pd[1659] = k44b;
pd[1660] = -Rap1Active*k44f;
pd[1663] = k43f;
pd[1670] = k20b + k21f;
pd[1681] = -MekInactive*k24f;
pd[1682] = -MekInactive*k26f;
pd[1696] = (-BRafActive*k26f - Raf1Active*k24f) + cj*(-1);
pd[1704] = k26b;
pd[1718] = k24b;
pd[1725] = k29f;
pd[1735] = MekInactive*k24f;
pd[1750] = Raf1Active*k24f;
pd[1772] = (-k24b - k25f) + cj*(-1);
pd[1786] = RasInactive*k10f;
pd[1801] = SosActive*k10f;
pd[1827] = (-k10b - k11f) + cj*(-1);
pd[1840] = -P90RskActive*k06f;
pd[1841] = (-SosActive*k06f) + cj*(-1);
pd[1868] = k06b + k07f;
pd[1872] = k09f;
pd[1900] = PP2AActive*k32f;
pd[1924] = ErkActive*k32f;
pd[1942] = (-k32b - k33f) + cj*(-1);
pd[1946] = PI3KInactive*k34f;
pd[1968] = boundEGFReceptor*k34f;
pd[1990] = (-k34b - k35f) + cj*(-1);
pd[2009] = (-AktInactive*k38f) + cj*(-1);
pd[2023] = -PI3KActive*k38f;
pd[2026] = k37f;
pd[2044] = k35f;
pd[2051] = k38b + k39f;
pd[2058] = -RasGapActive*k12f;
pd[2087] = k12b + k13f;
pd[2099] = (-RasActive*k12f) + cj*(-1);
pd[2114] = (-MekInactive*k26f - Raf1PPtase*k22f) + cj*(-1);
pd[2128] = -BRafActive*k26f;
pd[2136] = k26b + k27f;
pd[2139] = k22b;
pd[2148] = -BRafActive*k22f;
pd[2156] = k21f;
pd[2162] = SosInactive*k02f;
pd[2177] = boundEGFReceptor*k02f;
pd[2208] = (-k02b - k03f) + cj*(-1);
pd[2227] = -Rap1Inactive*k42f;
pd[2241] = (-C3GActive*k42f) + cj*(-1);
pd[2253] = k45f;
pd[2257] = k42b;
pd[2278] = (-P90RskInactive*k08f - PP2AActive*k32f) + cj*(-1);
pd[2286] = -ErkActive*k08f;
pd[2302] = -ErkActive*k32f;
pd[2304] = k08b + k09f;
pd[2306] = k31f;
pd[2320] = k32b;
pd[2323] = -freeNGFReceptor*k01f;
pd[2325] = k01b;
pd[2338] = (-NGF*k01f) + cj*(-1);
pd[2379] = -C3GInactive*k40f;
pd[2402] = (-boundNGFReceptor*k40f) + cj*(-1);
pd[2407] = k40b;
pd[2436] = -Raf1Inactive*k14f;
pd[2450] = (-RasActive*k14f) + cj*(-1);
pd[2459] = k17f;
pd[2467] = k14b;
pd[2471] = k19f;
pd[2486] = -PI3KInactive*k34f;
pd[2490] = -PI3KInactive*k36f;
pd[2508] = (-RasActive*k36f - boundEGFReceptor*k34f) + cj*(-1);
pd[2512] = k36b;
pd[2530] = k34b;
pd[2552] = -BRafInactive*k20f;
pd[2559] = (-Rap1Active*k20f) + cj*(-1);
pd[2571] = k23f;
pd[2588] = k20b;
pd[2606] = BRafInactive*k20f;
pd[2613] = Rap1Active*k20f;
pd[2642] = (-k20b - k21f) + cj*(-1);
pd[2652] = (-PI3KInactive*k36f - Raf1Inactive*k14f - RasGapActive*k12f) + cj*(-1);
pd[2666] = -RasActive*k14f;
pd[2670] = -RasActive*k36f;
pd[2674] = k36b + k37f;
pd[2681] = k12b;
pd[2683] = k14b + k15f;
pd[2691] = k11f;
pd[2693] = -RasActive*k12f;
pd[2704] = -RasInactive*k10f;
pd[2719] = (-SosActive*k10f) + cj*(-1);
pd[2735] = k13f;
pd[2745] = k10b;
pd[2763] = PP2AActive*k28f;
pd[2788] = MekActive*k28f;
pd[2805] = (-k28b - k29f) + cj*(-1);
pd[2819] = AktInactive*k38f;
pd[2833] = PI3KActive*k38f;
pd[2861] = (-k38b - k39f) + cj*(-1);
pd[2875] = (-Rap1Inactive*k42f) + cj*(-1);
pd[2889] = -C3GActive*k42f;
pd[2893] = k41f;
pd[2905] = k42b + k43f;
}
void ic_function(double *dynamicVars, double *constants) {
double cell = constants[0];
double RasGapTotal = constants[1];
double RapGapTotal = constants[2];
double PP2ATotal = constants[3];
double Raf1PPtaseTotal = constants[4];
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
double k00f = constants[21];
double k00b = constants[22];
double k01f = constants[23];
double k01b = constants[24];
double k02f = constants[25];
double k02b = constants[26];
double k03f = constants[27];
double k04f = constants[28];
double k04b = constants[29];
double k05f = constants[30];
double k06f = constants[31];
double k06b = constants[32];
double k07f = constants[33];
double k08f = constants[34];
double k08b = constants[35];
double k09f = constants[36];
double k10f = constants[37];
double k10b = constants[38];
double k11f = constants[39];
double k12f = constants[40];
double k12b = constants[41];
double k13f = constants[42];
double k14f = constants[43];
double k14b = constants[44];
double k15f = constants[45];
double k16f = constants[46];
double k16b = constants[47];
double k17f = constants[48];
double k18f = constants[49];
double k18b = constants[50];
double k19f = constants[51];
double k20f = constants[52];
double k20b = constants[53];
double k21f = constants[54];
double k22f = constants[55];
double k22b = constants[56];
double k23f = constants[57];
double k24f = constants[58];
double k24b = constants[59];
double k25f = constants[60];
double k26f = constants[61];
double k26b = constants[62];
double k27f = constants[63];
double k28f = constants[64];
double k28b = constants[65];
double k29f = constants[66];
double k30f = constants[67];
double k30b = constants[68];
double k31f = constants[69];
double k32f = constants[70];
double k32b = constants[71];
double k33f = constants[72];
double k34f = constants[73];
double k34b = constants[74];
double k35f = constants[75];
double k36f = constants[76];
double k36b = constants[77];
double k37f = constants[78];
double k38f = constants[79];
double k38b = constants[80];
double k39f = constants[81];
double k40f = constants[82];
double k40b = constants[83];
double k41f = constants[84];
double k42f = constants[85];
double k42b = constants[86];
double k43f = constants[87];
double k44f = constants[88];
double k44b = constants[89];
double k45f = constants[90];

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
dynamicVars[15] = Total_EGFReceptor;
dynamicVars[16] = Total_NGFReceptor;
dynamicVars[17] = Total_Sos;
dynamicVars[18] = Total_P90Rsk;
dynamicVars[19] = Total_Ras;
dynamicVars[20] = Total_Raf1;
dynamicVars[21] = Total_BRaf;
dynamicVars[22] = Total_Mek;
dynamicVars[23] = Total_Erk;
dynamicVars[24] = Total_PI3K;
dynamicVars[25] = Total_Akt;
dynamicVars[26] = Total_C3G;
dynamicVars[27] = Total_Rap1;
dynamicVars[28] = Zero;
dynamicVars[29] = Zero;
dynamicVars[30] = Zero;
dynamicVars[31] = Zero;
dynamicVars[32] = Zero;
dynamicVars[33] = Zero;
dynamicVars[34] = PP2ATotal;
dynamicVars[35] = Zero;
dynamicVars[36] = Zero;
dynamicVars[37] = Zero;
dynamicVars[38] = Zero;
dynamicVars[39] = Zero;
dynamicVars[40] = RapGapTotal;
dynamicVars[41] = Zero;
dynamicVars[42] = Raf1PPtaseTotal;
dynamicVars[43] = Zero;
dynamicVars[44] = Zero;
dynamicVars[45] = Zero;
dynamicVars[46] = Zero;
dynamicVars[47] = RasGapTotal;
dynamicVars[48] = Zero;
dynamicVars[49] = Zero;
dynamicVars[50] = Zero;
dynamicVars[51] = Zero;
dynamicVars[52] = Zero;
dynamicVars[53] = Zero;
}
