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
double dcell = constants[69];
double dRasGapActive = constants[70];
double dRapGapActive = constants[71];
double dPP2AActive = constants[72];
double dRaf1PPtase = constants[73];
double dEGF_IC = constants[74];
double dNGF_IC = constants[75];
double dTotal_EGFReceptor = constants[76];
double dTotal_NGFReceptor = constants[77];
double dTotal_Sos = constants[78];
double dTotal_P90Rsk = constants[79];
double dTotal_Ras = constants[80];
double dTotal_Raf1 = constants[81];
double dTotal_BRaf = constants[82];
double dTotal_Mek = constants[83];
double dTotal_Erk = constants[84];
double dTotal_PI3K = constants[85];
double dTotal_Akt = constants[86];
double dTotal_C3G = constants[87];
double dTotal_Rap1 = constants[88];
double dZero = constants[89];
double dkrbEGF = constants[90];
double dkruEGF = constants[91];
double dkrbNGF = constants[92];
double dkruNGF = constants[93];
double dkEGF = constants[94];
double dKmEGF = constants[95];
double dkNGF = constants[96];
double dKmNGF = constants[97];
double dkdSos = constants[98];
double dKmdSos = constants[99];
double dkSos = constants[100];
double dKmSos = constants[101];
double dkRasGap = constants[102];
double dKmRasGap = constants[103];
double dkRasToRaf1 = constants[104];
double dKmRasToRaf1 = constants[105];
double dkpRaf1 = constants[106];
double dKmpRaf1 = constants[107];
double dkpBRaf = constants[108];
double dKmpBRaf = constants[109];
double dkdMek = constants[110];
double dKmdMek = constants[111];
double dkpMekCytoplasmic = constants[112];
double dKmpMekCytoplasmic = constants[113];
double dkdErk = constants[114];
double dKmdErk = constants[115];
double dkpP90Rsk = constants[116];
double dKmpP90Rsk = constants[117];
double dkPI3K = constants[118];
double dKmPI3K = constants[119];
double dkPI3KRas = constants[120];
double dKmPI3KRas = constants[121];
double dkAkt = constants[122];
double dKmAkt = constants[123];
double dkdRaf1ByAkt = constants[124];
double dKmRaf1ByAkt = constants[125];
double dkC3GNGF = constants[126];
double dKmC3GNGF = constants[127];
double dkC3G = constants[128];
double dKmC3G = constants[129];
double dkRapGap = constants[130];
double dKmRapGap = constants[131];
double dkRap1ToBRaf = constants[132];
double dKmRap1ToBRaf = constants[133];
double dkdRaf1 = constants[134];
double dKmdRaf1 = constants[135];
double dkdBRaf = constants[136];
double dKmdBRaf = constants[137];

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
double dEGF = dynamicVars[15];
double dNGF = dynamicVars[16];
double dboundEGFReceptor = dynamicVars[17];
double dboundNGFReceptor = dynamicVars[18];
double dSosActive = dynamicVars[19];
double dP90RskActive = dynamicVars[20];
double dRasActive = dynamicVars[21];
double dRaf1Active = dynamicVars[22];
double dBRafActive = dynamicVars[23];
double dMekActive = dynamicVars[24];
double dErkActive = dynamicVars[25];
double dPI3KActive = dynamicVars[26];
double dAktActive = dynamicVars[27];
double dC3GActive = dynamicVars[28];
double dRap1Active = dynamicVars[29];
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
double dEGF_prime = yprime[15];
double dNGF_prime = yprime[16];
double dboundEGFReceptor_prime = yprime[17];
double dboundNGFReceptor_prime = yprime[18];
double dSosActive_prime = yprime[19];
double dP90RskActive_prime = yprime[20];
double dRasActive_prime = yprime[21];
double dRaf1Active_prime = yprime[22];
double dBRafActive_prime = yprime[23];
double dMekActive_prime = yprime[24];
double dErkActive_prime = yprime[25];
double dPI3KActive_prime = yprime[26];
double dAktActive_prime = yprime[27];
double dC3GActive_prime = yprime[28];
double dRap1Active_prime = yprime[29];

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
errors[15] = -EGF*cell*dkrbEGF*(Total_EGFReceptor - boundEGFReceptor) - EGF*cell*krbEGF*(dTotal_EGFReceptor - dboundEGFReceptor) - EGF*dcell*krbEGF*(Total_EGFReceptor - boundEGFReceptor) + boundEGFReceptor*cell*dkruEGF + boundEGFReceptor*dcell*kruEGF - cell*dEGF*krbEGF*(Total_EGFReceptor - boundEGFReceptor) + cell*dboundEGFReceptor*kruEGF - dEGF_prime;
errors[16] = -NGF*cell*dkrbNGF*(Total_NGFReceptor - boundNGFReceptor) - NGF*cell*krbNGF*(dTotal_NGFReceptor - dboundNGFReceptor) - NGF*dcell*krbNGF*(Total_NGFReceptor - boundNGFReceptor) + boundNGFReceptor*cell*dkruNGF + boundNGFReceptor*dcell*kruNGF - cell*dNGF*krbNGF*(Total_NGFReceptor - boundNGFReceptor) + cell*dboundNGFReceptor*kruNGF - dNGF_prime;
errors[17] = EGF*cell*dkrbEGF*(Total_EGFReceptor - boundEGFReceptor) + EGF*cell*krbEGF*(dTotal_EGFReceptor - dboundEGFReceptor) + EGF*dcell*krbEGF*(Total_EGFReceptor - boundEGFReceptor) - boundEGFReceptor*cell*dkruEGF - boundEGFReceptor*dcell*kruEGF + cell*dEGF*krbEGF*(Total_EGFReceptor - boundEGFReceptor) - cell*dboundEGFReceptor*kruEGF - dboundEGFReceptor_prime;
errors[18] = NGF*cell*dkrbNGF*(Total_NGFReceptor - boundNGFReceptor) + NGF*cell*krbNGF*(dTotal_NGFReceptor - dboundNGFReceptor) + NGF*dcell*krbNGF*(Total_NGFReceptor - boundNGFReceptor) - boundNGFReceptor*cell*dkruNGF - boundNGFReceptor*dcell*kruNGF + cell*dNGF*krbNGF*(Total_NGFReceptor - boundNGFReceptor) - cell*dboundNGFReceptor*kruNGF - dboundNGFReceptor_prime;
errors[19] = -P90RskActive*SosActive*cell*dkdSos/(KmdSos + SosActive) - P90RskActive*SosActive*cell*kdSos*(-dKmdSos - dSosActive)/pow(KmdSos + SosActive, 2) - P90RskActive*SosActive*dcell*kdSos/(KmdSos + SosActive) - P90RskActive*cell*dSosActive*kdSos/(KmdSos + SosActive) - SosActive*cell*dP90RskActive*kdSos/(KmdSos + SosActive) + boundEGFReceptor*cell*dkEGF*(-SosActive + Total_Sos)/(KmEGF - SosActive + Total_Sos) + boundEGFReceptor*cell*kEGF*(-SosActive + Total_Sos)*(-dKmEGF + dSosActive - dTotal_Sos)/pow(KmEGF - SosActive + Total_Sos, 2) + boundEGFReceptor*cell*kEGF*(-dSosActive + dTotal_Sos)/(KmEGF - SosActive + Total_Sos) + boundEGFReceptor*dcell*kEGF*(-SosActive + Total_Sos)/(KmEGF - SosActive + Total_Sos) + boundNGFReceptor*cell*dkNGF*(-SosActive + Total_Sos)/(KmNGF - SosActive + Total_Sos) + boundNGFReceptor*cell*kNGF*(-SosActive + Total_Sos)*(-dKmNGF + dSosActive - dTotal_Sos)/pow(KmNGF - SosActive + Total_Sos, 2) + boundNGFReceptor*cell*kNGF*(-dSosActive + dTotal_Sos)/(KmNGF - SosActive + Total_Sos) + boundNGFReceptor*dcell*kNGF*(-SosActive + Total_Sos)/(KmNGF - SosActive + Total_Sos) + cell*dboundEGFReceptor*kEGF*(-SosActive + Total_Sos)/(KmEGF - SosActive + Total_Sos) + cell*dboundNGFReceptor*kNGF*(-SosActive + Total_Sos)/(KmNGF - SosActive + Total_Sos) - dSosActive_prime;
errors[20] = ErkActive*cell*dkpP90Rsk*(-P90RskActive + Total_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk) + ErkActive*cell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)*(-dKmpP90Rsk + dP90RskActive - dTotal_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 2) + ErkActive*cell*kpP90Rsk*(-dP90RskActive + dTotal_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk) + ErkActive*dcell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk) + cell*dErkActive*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk) - dP90RskActive_prime;
errors[21] = -RasActive*RasGapActive*cell*dkRasGap/(KmRasGap + RasActive) - RasActive*RasGapActive*cell*kRasGap*(-dKmRasGap - dRasActive)/pow(KmRasGap + RasActive, 2) - RasActive*RasGapActive*dcell*kRasGap/(KmRasGap + RasActive) - RasActive*cell*dRasGapActive*kRasGap/(KmRasGap + RasActive) - RasGapActive*cell*dRasActive*kRasGap/(KmRasGap + RasActive) + SosActive*cell*dkSos*(-RasActive + Total_Ras)/(KmSos - RasActive + Total_Ras) + SosActive*cell*kSos*(-RasActive + Total_Ras)*(-dKmSos + dRasActive - dTotal_Ras)/pow(KmSos - RasActive + Total_Ras, 2) + SosActive*cell*kSos*(-dRasActive + dTotal_Ras)/(KmSos - RasActive + Total_Ras) + SosActive*dcell*kSos*(-RasActive + Total_Ras)/(KmSos - RasActive + Total_Ras) + cell*dSosActive*kSos*(-RasActive + Total_Ras)/(KmSos - RasActive + Total_Ras) - dRasActive_prime;
errors[22] = -AktActive*Raf1Active*cell*dkdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) - AktActive*Raf1Active*cell*kdRaf1ByAkt*(-dKmRaf1ByAkt - dRaf1Active)/pow(KmRaf1ByAkt + Raf1Active, 2) - AktActive*Raf1Active*dcell*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) - AktActive*cell*dRaf1Active*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) - Raf1Active*Raf1PPtase*cell*dkdRaf1/(KmdRaf1 + Raf1Active) - Raf1Active*Raf1PPtase*cell*kdRaf1*(-dKmdRaf1 - dRaf1Active)/pow(KmdRaf1 + Raf1Active, 2) - Raf1Active*Raf1PPtase*dcell*kdRaf1/(KmdRaf1 + Raf1Active) - Raf1Active*cell*dAktActive*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) - Raf1Active*cell*dRaf1PPtase*kdRaf1/(KmdRaf1 + Raf1Active) - Raf1PPtase*cell*dRaf1Active*kdRaf1/(KmdRaf1 + Raf1Active) + RasActive*cell*dkRasToRaf1*(-Raf1Active + Total_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1) + RasActive*cell*kRasToRaf1*(-Raf1Active + Total_Raf1)*(-dKmRasToRaf1 + dRaf1Active - dTotal_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 2) + RasActive*cell*kRasToRaf1*(-dRaf1Active + dTotal_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1) + RasActive*dcell*kRasToRaf1*(-Raf1Active + Total_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1) + cell*dRasActive*kRasToRaf1*(-Raf1Active + Total_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1) - dRaf1Active_prime;
errors[23] = -BRafActive*Raf1PPtase*cell*dkdBRaf/(BRafActive + KmdBRaf) - BRafActive*Raf1PPtase*cell*kdBRaf*(-dBRafActive - dKmdBRaf)/pow(BRafActive + KmdBRaf, 2) - BRafActive*Raf1PPtase*dcell*kdBRaf/(BRafActive + KmdBRaf) - BRafActive*cell*dRaf1PPtase*kdBRaf/(BRafActive + KmdBRaf) - Raf1PPtase*cell*dBRafActive*kdBRaf/(BRafActive + KmdBRaf) + Rap1Active*cell*dkRap1ToBRaf*(-BRafActive + Total_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf) + Rap1Active*cell*kRap1ToBRaf*(-BRafActive + Total_BRaf)*(dBRafActive - dKmRap1ToBRaf - dTotal_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 2) + Rap1Active*cell*kRap1ToBRaf*(-dBRafActive + dTotal_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf) + Rap1Active*dcell*kRap1ToBRaf*(-BRafActive + Total_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf) + cell*dRap1Active*kRap1ToBRaf*(-BRafActive + Total_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf) - dBRafActive_prime;
errors[24] = BRafActive*cell*dkpBRaf*(-MekActive + Total_Mek)/(KmpBRaf - MekActive + Total_Mek) + BRafActive*cell*kpBRaf*(-MekActive + Total_Mek)*(-dKmpBRaf + dMekActive - dTotal_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 2) + BRafActive*cell*kpBRaf*(-dMekActive + dTotal_Mek)/(KmpBRaf - MekActive + Total_Mek) + BRafActive*dcell*kpBRaf*(-MekActive + Total_Mek)/(KmpBRaf - MekActive + Total_Mek) - MekActive*PP2AActive*cell*dkdMek/(KmdMek + MekActive) - MekActive*PP2AActive*cell*kdMek*(-dKmdMek - dMekActive)/pow(KmdMek + MekActive, 2) - MekActive*PP2AActive*dcell*kdMek/(KmdMek + MekActive) - MekActive*cell*dPP2AActive*kdMek/(KmdMek + MekActive) - PP2AActive*cell*dMekActive*kdMek/(KmdMek + MekActive) + Raf1Active*cell*dkpRaf1*(-MekActive + Total_Mek)/(KmpRaf1 - MekActive + Total_Mek) + Raf1Active*cell*kpRaf1*(-MekActive + Total_Mek)*(-dKmpRaf1 + dMekActive - dTotal_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 2) + Raf1Active*cell*kpRaf1*(-dMekActive + dTotal_Mek)/(KmpRaf1 - MekActive + Total_Mek) + Raf1Active*dcell*kpRaf1*(-MekActive + Total_Mek)/(KmpRaf1 - MekActive + Total_Mek) + cell*dBRafActive*kpBRaf*(-MekActive + Total_Mek)/(KmpBRaf - MekActive + Total_Mek) + cell*dRaf1Active*kpRaf1*(-MekActive + Total_Mek)/(KmpRaf1 - MekActive + Total_Mek) - dMekActive_prime;
errors[25] = -ErkActive*PP2AActive*cell*dkdErk/(ErkActive + KmdErk) - ErkActive*PP2AActive*cell*kdErk*(-dErkActive - dKmdErk)/pow(ErkActive + KmdErk, 2) - ErkActive*PP2AActive*dcell*kdErk/(ErkActive + KmdErk) - ErkActive*cell*dPP2AActive*kdErk/(ErkActive + KmdErk) + MekActive*cell*dkpMekCytoplasmic*(-ErkActive + Total_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk) + MekActive*cell*kpMekCytoplasmic*(-ErkActive + Total_Erk)*(dErkActive - dKmpMekCytoplasmic - dTotal_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 2) + MekActive*cell*kpMekCytoplasmic*(-dErkActive + dTotal_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk) + MekActive*dcell*kpMekCytoplasmic*(-ErkActive + Total_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk) - PP2AActive*cell*dErkActive*kdErk/(ErkActive + KmdErk) + cell*dMekActive*kpMekCytoplasmic*(-ErkActive + Total_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk) - dErkActive_prime;
errors[26] = RasActive*cell*dkPI3KRas*(-PI3KActive + Total_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K) + RasActive*cell*kPI3KRas*(-PI3KActive + Total_PI3K)*(-dKmPI3KRas + dPI3KActive - dTotal_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 2) + RasActive*cell*kPI3KRas*(-dPI3KActive + dTotal_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K) + RasActive*dcell*kPI3KRas*(-PI3KActive + Total_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K) + boundEGFReceptor*cell*dkPI3K*(-PI3KActive + Total_PI3K)/(KmPI3K - PI3KActive + Total_PI3K) + boundEGFReceptor*cell*kPI3K*(-PI3KActive + Total_PI3K)*(-dKmPI3K + dPI3KActive - dTotal_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 2) + boundEGFReceptor*cell*kPI3K*(-dPI3KActive + dTotal_PI3K)/(KmPI3K - PI3KActive + Total_PI3K) + boundEGFReceptor*dcell*kPI3K*(-PI3KActive + Total_PI3K)/(KmPI3K - PI3KActive + Total_PI3K) + cell*dRasActive*kPI3KRas*(-PI3KActive + Total_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K) + cell*dboundEGFReceptor*kPI3K*(-PI3KActive + Total_PI3K)/(KmPI3K - PI3KActive + Total_PI3K) - dPI3KActive_prime;
errors[27] = PI3KActive*cell*dkAkt*(-AktActive + Total_Akt)/(-AktActive + KmAkt + Total_Akt) + PI3KActive*cell*kAkt*(-AktActive + Total_Akt)*(dAktActive - dKmAkt - dTotal_Akt)/pow(-AktActive + KmAkt + Total_Akt, 2) + PI3KActive*cell*kAkt*(-dAktActive + dTotal_Akt)/(-AktActive + KmAkt + Total_Akt) + PI3KActive*dcell*kAkt*(-AktActive + Total_Akt)/(-AktActive + KmAkt + Total_Akt) + cell*dPI3KActive*kAkt*(-AktActive + Total_Akt)/(-AktActive + KmAkt + Total_Akt) - dAktActive_prime;
errors[28] = boundNGFReceptor*cell*dkC3GNGF*(-C3GActive + Total_C3G)/(-C3GActive + KmC3GNGF + Total_C3G) + boundNGFReceptor*cell*kC3GNGF*(-C3GActive + Total_C3G)*(dC3GActive - dKmC3GNGF - dTotal_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 2) + boundNGFReceptor*cell*kC3GNGF*(-dC3GActive + dTotal_C3G)/(-C3GActive + KmC3GNGF + Total_C3G) + boundNGFReceptor*dcell*kC3GNGF*(-C3GActive + Total_C3G)/(-C3GActive + KmC3GNGF + Total_C3G) + cell*dboundNGFReceptor*kC3GNGF*(-C3GActive + Total_C3G)/(-C3GActive + KmC3GNGF + Total_C3G) - dC3GActive_prime;
errors[29] = C3GActive*cell*dkC3G*(-Rap1Active + Total_Rap1)/(KmC3G - Rap1Active + Total_Rap1) + C3GActive*cell*kC3G*(-Rap1Active + Total_Rap1)*(-dKmC3G + dRap1Active - dTotal_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 2) + C3GActive*cell*kC3G*(-dRap1Active + dTotal_Rap1)/(KmC3G - Rap1Active + Total_Rap1) + C3GActive*dcell*kC3G*(-Rap1Active + Total_Rap1)/(KmC3G - Rap1Active + Total_Rap1) - Rap1Active*RapGapActive*cell*dkRapGap/(KmRapGap + Rap1Active) - Rap1Active*RapGapActive*cell*kRapGap*(-dKmRapGap - dRap1Active)/pow(KmRapGap + Rap1Active, 2) - Rap1Active*RapGapActive*dcell*kRapGap/(KmRapGap + Rap1Active) - Rap1Active*cell*dRapGapActive*kRapGap/(KmRapGap + Rap1Active) - RapGapActive*cell*dRap1Active*kRapGap/(KmRapGap + Rap1Active) + cell*dC3GActive*kC3G*(-Rap1Active + Total_Rap1)/(KmC3G - Rap1Active + Total_Rap1) - dRap1Active_prime;
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
double dcell = constants[69];
double dRasGapActive = constants[70];
double dRapGapActive = constants[71];
double dPP2AActive = constants[72];
double dRaf1PPtase = constants[73];
double dEGF_IC = constants[74];
double dNGF_IC = constants[75];
double dTotal_EGFReceptor = constants[76];
double dTotal_NGFReceptor = constants[77];
double dTotal_Sos = constants[78];
double dTotal_P90Rsk = constants[79];
double dTotal_Ras = constants[80];
double dTotal_Raf1 = constants[81];
double dTotal_BRaf = constants[82];
double dTotal_Mek = constants[83];
double dTotal_Erk = constants[84];
double dTotal_PI3K = constants[85];
double dTotal_Akt = constants[86];
double dTotal_C3G = constants[87];
double dTotal_Rap1 = constants[88];
double dZero = constants[89];
double dkrbEGF = constants[90];
double dkruEGF = constants[91];
double dkrbNGF = constants[92];
double dkruNGF = constants[93];
double dkEGF = constants[94];
double dKmEGF = constants[95];
double dkNGF = constants[96];
double dKmNGF = constants[97];
double dkdSos = constants[98];
double dKmdSos = constants[99];
double dkSos = constants[100];
double dKmSos = constants[101];
double dkRasGap = constants[102];
double dKmRasGap = constants[103];
double dkRasToRaf1 = constants[104];
double dKmRasToRaf1 = constants[105];
double dkpRaf1 = constants[106];
double dKmpRaf1 = constants[107];
double dkpBRaf = constants[108];
double dKmpBRaf = constants[109];
double dkdMek = constants[110];
double dKmdMek = constants[111];
double dkpMekCytoplasmic = constants[112];
double dKmpMekCytoplasmic = constants[113];
double dkdErk = constants[114];
double dKmdErk = constants[115];
double dkpP90Rsk = constants[116];
double dKmpP90Rsk = constants[117];
double dkPI3K = constants[118];
double dKmPI3K = constants[119];
double dkPI3KRas = constants[120];
double dKmPI3KRas = constants[121];
double dkAkt = constants[122];
double dKmAkt = constants[123];
double dkdRaf1ByAkt = constants[124];
double dKmRaf1ByAkt = constants[125];
double dkC3GNGF = constants[126];
double dKmC3GNGF = constants[127];
double dkC3G = constants[128];
double dKmC3G = constants[129];
double dkRapGap = constants[130];
double dKmRapGap = constants[131];
double dkRap1ToBRaf = constants[132];
double dKmRap1ToBRaf = constants[133];
double dkdRaf1 = constants[134];
double dKmdRaf1 = constants[135];
double dkdBRaf = constants[136];
double dKmdBRaf = constants[137];

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
double dEGF = dynamicVars[15];
double dNGF = dynamicVars[16];
double dboundEGFReceptor = dynamicVars[17];
double dboundNGFReceptor = dynamicVars[18];
double dSosActive = dynamicVars[19];
double dP90RskActive = dynamicVars[20];
double dRasActive = dynamicVars[21];
double dRaf1Active = dynamicVars[22];
double dBRafActive = dynamicVars[23];
double dMekActive = dynamicVars[24];
double dErkActive = dynamicVars[25];
double dPI3KActive = dynamicVars[26];
double dAktActive = dynamicVars[27];
double dC3GActive = dynamicVars[28];
double dRap1Active = dynamicVars[29];

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
double dEGF_prime = yprime[15];
double dNGF_prime = yprime[16];
double dboundEGFReceptor_prime = yprime[17];
double dboundNGFReceptor_prime = yprime[18];
double dSosActive_prime = yprime[19];
double dP90RskActive_prime = yprime[20];
double dRasActive_prime = yprime[21];
double dRaf1Active_prime = yprime[22];
double dBRafActive_prime = yprime[23];
double dMekActive_prime = yprime[24];
double dErkActive_prime = yprime[25];
double dPI3KActive_prime = yprime[26];
double dAktActive_prime = yprime[27];
double dC3GActive_prime = yprime[28];
double dRap1Active_prime = yprime[29];

memset(pd, 0, sizeof(pd)*900);
pd[0] = (-cell*krbEGF*(Total_EGFReceptor - boundEGFReceptor)) + cj*(-1);
pd[2] = EGF*cell*krbEGF + cell*kruEGF;
pd[31] = (-cell*krbNGF*(Total_NGFReceptor - boundNGFReceptor)) + cj*(-1);
pd[33] = NGF*cell*krbNGF + cell*kruNGF;
pd[60] = cell*krbEGF*(Total_EGFReceptor - boundEGFReceptor);
pd[62] = (-EGF*cell*krbEGF - cell*kruEGF) + cj*(-1);
pd[91] = cell*krbNGF*(Total_NGFReceptor - boundNGFReceptor);
pd[93] = (-NGF*cell*krbNGF - cell*kruNGF) + cj*(-1);
pd[122] = cell*kEGF*(-SosActive + Total_Sos)/(KmEGF - SosActive + Total_Sos);
pd[123] = cell*kNGF*(-SosActive + Total_Sos)/(KmNGF - SosActive + Total_Sos);
pd[124] = (P90RskActive*SosActive*cell*kdSos/pow(KmdSos + SosActive, 2) - P90RskActive*cell*kdSos/(KmdSos + SosActive) + boundEGFReceptor*cell*kEGF*(-SosActive + Total_Sos)/pow(KmEGF - SosActive + Total_Sos, 2) - boundEGFReceptor*cell*kEGF/(KmEGF - SosActive + Total_Sos) + boundNGFReceptor*cell*kNGF*(-SosActive + Total_Sos)/pow(KmNGF - SosActive + Total_Sos, 2) - boundNGFReceptor*cell*kNGF/(KmNGF - SosActive + Total_Sos)) + cj*(-1);
pd[125] = -SosActive*cell*kdSos/(KmdSos + SosActive);
pd[155] = (ErkActive*cell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 2) - ErkActive*cell*kpP90Rsk/(KmpP90Rsk - P90RskActive + Total_P90Rsk)) + cj*(-1);
pd[160] = cell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk);
pd[184] = cell*kSos*(-RasActive + Total_Ras)/(KmSos - RasActive + Total_Ras);
pd[186] = (RasActive*RasGapActive*cell*kRasGap/pow(KmRasGap + RasActive, 2) - RasGapActive*cell*kRasGap/(KmRasGap + RasActive) + SosActive*cell*kSos*(-RasActive + Total_Ras)/pow(KmSos - RasActive + Total_Ras, 2) - SosActive*cell*kSos/(KmSos - RasActive + Total_Ras)) + cj*(-1);
pd[216] = cell*kRasToRaf1*(-Raf1Active + Total_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1);
pd[217] = (AktActive*Raf1Active*cell*kdRaf1ByAkt/pow(KmRaf1ByAkt + Raf1Active, 2) - AktActive*cell*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) + Raf1Active*Raf1PPtase*cell*kdRaf1/pow(KmdRaf1 + Raf1Active, 2) - Raf1PPtase*cell*kdRaf1/(KmdRaf1 + Raf1Active) + RasActive*cell*kRasToRaf1*(-Raf1Active + Total_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 2) - RasActive*cell*kRasToRaf1/(KmRasToRaf1 - Raf1Active + Total_Raf1)) + cj*(-1);
pd[222] = -Raf1Active*cell*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active);
pd[248] = (BRafActive*Raf1PPtase*cell*kdBRaf/pow(BRafActive + KmdBRaf, 2) - Raf1PPtase*cell*kdBRaf/(BRafActive + KmdBRaf) + Rap1Active*cell*kRap1ToBRaf*(-BRafActive + Total_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 2) - Rap1Active*cell*kRap1ToBRaf/(-BRafActive + KmRap1ToBRaf + Total_BRaf)) + cj*(-1);
pd[254] = cell*kRap1ToBRaf*(-BRafActive + Total_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf);
pd[277] = cell*kpRaf1*(-MekActive + Total_Mek)/(KmpRaf1 - MekActive + Total_Mek);
pd[278] = cell*kpBRaf*(-MekActive + Total_Mek)/(KmpBRaf - MekActive + Total_Mek);
pd[279] = (BRafActive*cell*kpBRaf*(-MekActive + Total_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 2) - BRafActive*cell*kpBRaf/(KmpBRaf - MekActive + Total_Mek) + MekActive*PP2AActive*cell*kdMek/pow(KmdMek + MekActive, 2) - PP2AActive*cell*kdMek/(KmdMek + MekActive) + Raf1Active*cell*kpRaf1*(-MekActive + Total_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 2) - Raf1Active*cell*kpRaf1/(KmpRaf1 - MekActive + Total_Mek)) + cj*(-1);
pd[309] = cell*kpMekCytoplasmic*(-ErkActive + Total_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk);
pd[310] = (ErkActive*PP2AActive*cell*kdErk/pow(ErkActive + KmdErk, 2) + MekActive*cell*kpMekCytoplasmic*(-ErkActive + Total_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 2) - MekActive*cell*kpMekCytoplasmic/(-ErkActive + KmpMekCytoplasmic + Total_Erk) - PP2AActive*cell*kdErk/(ErkActive + KmdErk)) + cj*(-1);
pd[332] = cell*kPI3K*(-PI3KActive + Total_PI3K)/(KmPI3K - PI3KActive + Total_PI3K);
pd[336] = cell*kPI3KRas*(-PI3KActive + Total_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K);
pd[341] = (RasActive*cell*kPI3KRas*(-PI3KActive + Total_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 2) - RasActive*cell*kPI3KRas/(KmPI3KRas - PI3KActive + Total_PI3K) + boundEGFReceptor*cell*kPI3K*(-PI3KActive + Total_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 2) - boundEGFReceptor*cell*kPI3K/(KmPI3K - PI3KActive + Total_PI3K)) + cj*(-1);
pd[371] = cell*kAkt*(-AktActive + Total_Akt)/(-AktActive + KmAkt + Total_Akt);
pd[372] = (PI3KActive*cell*kAkt*(-AktActive + Total_Akt)/pow(-AktActive + KmAkt + Total_Akt, 2) - PI3KActive*cell*kAkt/(-AktActive + KmAkt + Total_Akt)) + cj*(-1);
pd[393] = cell*kC3GNGF*(-C3GActive + Total_C3G)/(-C3GActive + KmC3GNGF + Total_C3G);
pd[403] = (boundNGFReceptor*cell*kC3GNGF*(-C3GActive + Total_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 2) - boundNGFReceptor*cell*kC3GNGF/(-C3GActive + KmC3GNGF + Total_C3G)) + cj*(-1);
pd[433] = cell*kC3G*(-Rap1Active + Total_Rap1)/(KmC3G - Rap1Active + Total_Rap1);
pd[434] = (C3GActive*cell*kC3G*(-Rap1Active + Total_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 2) - C3GActive*cell*kC3G/(KmC3G - Rap1Active + Total_Rap1) + Rap1Active*RapGapActive*cell*kRapGap/pow(KmRapGap + Rap1Active, 2) - RapGapActive*cell*kRapGap/(KmRapGap + Rap1Active)) + cj*(-1);
pd[450] = -cell*dkrbEGF*(Total_EGFReceptor - boundEGFReceptor) - cell*krbEGF*(dTotal_EGFReceptor - dboundEGFReceptor) - dcell*krbEGF*(Total_EGFReceptor - boundEGFReceptor);
pd[452] = EGF*cell*dkrbEGF + EGF*dcell*krbEGF + cell*dEGF*krbEGF + cell*dkruEGF + dcell*kruEGF;
pd[465] = (-cell*krbEGF*(Total_EGFReceptor - boundEGFReceptor)) + cj*(-1);
pd[467] = EGF*cell*krbEGF + cell*kruEGF;
pd[481] = -cell*dkrbNGF*(Total_NGFReceptor - boundNGFReceptor) - cell*krbNGF*(dTotal_NGFReceptor - dboundNGFReceptor) - dcell*krbNGF*(Total_NGFReceptor - boundNGFReceptor);
pd[483] = NGF*cell*dkrbNGF + NGF*dcell*krbNGF + cell*dNGF*krbNGF + cell*dkruNGF + dcell*kruNGF;
pd[496] = (-cell*krbNGF*(Total_NGFReceptor - boundNGFReceptor)) + cj*(-1);
pd[498] = NGF*cell*krbNGF + cell*kruNGF;
pd[510] = cell*dkrbEGF*(Total_EGFReceptor - boundEGFReceptor) + cell*krbEGF*(dTotal_EGFReceptor - dboundEGFReceptor) + dcell*krbEGF*(Total_EGFReceptor - boundEGFReceptor);
pd[512] = -EGF*cell*dkrbEGF - EGF*dcell*krbEGF - cell*dEGF*krbEGF - cell*dkruEGF - dcell*kruEGF;
pd[525] = cell*krbEGF*(Total_EGFReceptor - boundEGFReceptor);
pd[527] = (-EGF*cell*krbEGF - cell*kruEGF) + cj*(-1);
pd[541] = cell*dkrbNGF*(Total_NGFReceptor - boundNGFReceptor) + cell*krbNGF*(dTotal_NGFReceptor - dboundNGFReceptor) + dcell*krbNGF*(Total_NGFReceptor - boundNGFReceptor);
pd[543] = -NGF*cell*dkrbNGF - NGF*dcell*krbNGF - cell*dNGF*krbNGF - cell*dkruNGF - dcell*kruNGF;
pd[556] = cell*krbNGF*(Total_NGFReceptor - boundNGFReceptor);
pd[558] = (-NGF*cell*krbNGF - cell*kruNGF) + cj*(-1);
pd[572] = cell*dkEGF*(-SosActive + Total_Sos)/(KmEGF - SosActive + Total_Sos) + cell*kEGF*(-SosActive + Total_Sos)*(-dKmEGF + dSosActive - dTotal_Sos)/pow(KmEGF - SosActive + Total_Sos, 2) + cell*kEGF*(-dSosActive + dTotal_Sos)/(KmEGF - SosActive + Total_Sos) + dcell*kEGF*(-SosActive + Total_Sos)/(KmEGF - SosActive + Total_Sos);
pd[573] = cell*dkNGF*(-SosActive + Total_Sos)/(KmNGF - SosActive + Total_Sos) + cell*kNGF*(-SosActive + Total_Sos)*(-dKmNGF + dSosActive - dTotal_Sos)/pow(KmNGF - SosActive + Total_Sos, 2) + cell*kNGF*(-dSosActive + dTotal_Sos)/(KmNGF - SosActive + Total_Sos) + dcell*kNGF*(-SosActive + Total_Sos)/(KmNGF - SosActive + Total_Sos);
pd[574] = P90RskActive*SosActive*cell*dkdSos/pow(KmdSos + SosActive, 2) + 2*P90RskActive*SosActive*cell*kdSos*(-dKmdSos - dSosActive)/pow(KmdSos + SosActive, 3) + P90RskActive*SosActive*dcell*kdSos/pow(KmdSos + SosActive, 2) + P90RskActive*cell*dSosActive*kdSos/pow(KmdSos + SosActive, 2) - P90RskActive*cell*dkdSos/(KmdSos + SosActive) - P90RskActive*cell*kdSos*(-dKmdSos - dSosActive)/pow(KmdSos + SosActive, 2) - P90RskActive*dcell*kdSos/(KmdSos + SosActive) + SosActive*cell*dP90RskActive*kdSos/pow(KmdSos + SosActive, 2) + boundEGFReceptor*cell*dkEGF*(-SosActive + Total_Sos)/pow(KmEGF - SosActive + Total_Sos, 2) - boundEGFReceptor*cell*dkEGF/(KmEGF - SosActive + Total_Sos) + 2*boundEGFReceptor*cell*kEGF*(-SosActive + Total_Sos)*(-dKmEGF + dSosActive - dTotal_Sos)/pow(KmEGF - SosActive + Total_Sos, 3) + boundEGFReceptor*cell*kEGF*(-dSosActive + dTotal_Sos)/pow(KmEGF - SosActive + Total_Sos, 2) - boundEGFReceptor*cell*kEGF*(-dKmEGF + dSosActive - dTotal_Sos)/pow(KmEGF - SosActive + Total_Sos, 2) + boundEGFReceptor*dcell*kEGF*(-SosActive + Total_Sos)/pow(KmEGF - SosActive + Total_Sos, 2) - boundEGFReceptor*dcell*kEGF/(KmEGF - SosActive + Total_Sos) + boundNGFReceptor*cell*dkNGF*(-SosActive + Total_Sos)/pow(KmNGF - SosActive + Total_Sos, 2) - boundNGFReceptor*cell*dkNGF/(KmNGF - SosActive + Total_Sos) + 2*boundNGFReceptor*cell*kNGF*(-SosActive + Total_Sos)*(-dKmNGF + dSosActive - dTotal_Sos)/pow(KmNGF - SosActive + Total_Sos, 3) + boundNGFReceptor*cell*kNGF*(-dSosActive + dTotal_Sos)/pow(KmNGF - SosActive + Total_Sos, 2) - boundNGFReceptor*cell*kNGF*(-dKmNGF + dSosActive - dTotal_Sos)/pow(KmNGF - SosActive + Total_Sos, 2) + boundNGFReceptor*dcell*kNGF*(-SosActive + Total_Sos)/pow(KmNGF - SosActive + Total_Sos, 2) - boundNGFReceptor*dcell*kNGF/(KmNGF - SosActive + Total_Sos) - cell*dP90RskActive*kdSos/(KmdSos + SosActive) + cell*dboundEGFReceptor*kEGF*(-SosActive + Total_Sos)/pow(KmEGF - SosActive + Total_Sos, 2) - cell*dboundEGFReceptor*kEGF/(KmEGF - SosActive + Total_Sos) + cell*dboundNGFReceptor*kNGF*(-SosActive + Total_Sos)/pow(KmNGF - SosActive + Total_Sos, 2) - cell*dboundNGFReceptor*kNGF/(KmNGF - SosActive + Total_Sos);
pd[575] = -SosActive*cell*dkdSos/(KmdSos + SosActive) - SosActive*cell*kdSos*(-dKmdSos - dSosActive)/pow(KmdSos + SosActive, 2) - SosActive*dcell*kdSos/(KmdSos + SosActive) - cell*dSosActive*kdSos/(KmdSos + SosActive);
pd[587] = cell*kEGF*(-SosActive + Total_Sos)/(KmEGF - SosActive + Total_Sos);
pd[588] = cell*kNGF*(-SosActive + Total_Sos)/(KmNGF - SosActive + Total_Sos);
pd[589] = (P90RskActive*SosActive*cell*kdSos/pow(KmdSos + SosActive, 2) - P90RskActive*cell*kdSos/(KmdSos + SosActive) + boundEGFReceptor*cell*kEGF*(-SosActive + Total_Sos)/pow(KmEGF - SosActive + Total_Sos, 2) - boundEGFReceptor*cell*kEGF/(KmEGF - SosActive + Total_Sos) + boundNGFReceptor*cell*kNGF*(-SosActive + Total_Sos)/pow(KmNGF - SosActive + Total_Sos, 2) - boundNGFReceptor*cell*kNGF/(KmNGF - SosActive + Total_Sos)) + cj*(-1);
pd[590] = -SosActive*cell*kdSos/(KmdSos + SosActive);
pd[605] = ErkActive*cell*dkpP90Rsk*(-P90RskActive + Total_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 2) - ErkActive*cell*dkpP90Rsk/(KmpP90Rsk - P90RskActive + Total_P90Rsk) + 2*ErkActive*cell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)*(-dKmpP90Rsk + dP90RskActive - dTotal_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 3) + ErkActive*cell*kpP90Rsk*(-dP90RskActive + dTotal_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 2) - ErkActive*cell*kpP90Rsk*(-dKmpP90Rsk + dP90RskActive - dTotal_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 2) + ErkActive*dcell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 2) - ErkActive*dcell*kpP90Rsk/(KmpP90Rsk - P90RskActive + Total_P90Rsk) + cell*dErkActive*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 2) - cell*dErkActive*kpP90Rsk/(KmpP90Rsk - P90RskActive + Total_P90Rsk);
pd[610] = cell*dkpP90Rsk*(-P90RskActive + Total_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk) + cell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)*(-dKmpP90Rsk + dP90RskActive - dTotal_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 2) + cell*kpP90Rsk*(-dP90RskActive + dTotal_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk) + dcell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk);
pd[620] = (ErkActive*cell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/pow(KmpP90Rsk - P90RskActive + Total_P90Rsk, 2) - ErkActive*cell*kpP90Rsk/(KmpP90Rsk - P90RskActive + Total_P90Rsk)) + cj*(-1);
pd[625] = cell*kpP90Rsk*(-P90RskActive + Total_P90Rsk)/(KmpP90Rsk - P90RskActive + Total_P90Rsk);
pd[634] = cell*dkSos*(-RasActive + Total_Ras)/(KmSos - RasActive + Total_Ras) + cell*kSos*(-RasActive + Total_Ras)*(-dKmSos + dRasActive - dTotal_Ras)/pow(KmSos - RasActive + Total_Ras, 2) + cell*kSos*(-dRasActive + dTotal_Ras)/(KmSos - RasActive + Total_Ras) + dcell*kSos*(-RasActive + Total_Ras)/(KmSos - RasActive + Total_Ras);
pd[636] = RasActive*RasGapActive*cell*dkRasGap/pow(KmRasGap + RasActive, 2) + 2*RasActive*RasGapActive*cell*kRasGap*(-dKmRasGap - dRasActive)/pow(KmRasGap + RasActive, 3) + RasActive*RasGapActive*dcell*kRasGap/pow(KmRasGap + RasActive, 2) + RasActive*cell*dRasGapActive*kRasGap/pow(KmRasGap + RasActive, 2) + RasGapActive*cell*dRasActive*kRasGap/pow(KmRasGap + RasActive, 2) - RasGapActive*cell*dkRasGap/(KmRasGap + RasActive) - RasGapActive*cell*kRasGap*(-dKmRasGap - dRasActive)/pow(KmRasGap + RasActive, 2) - RasGapActive*dcell*kRasGap/(KmRasGap + RasActive) + SosActive*cell*dkSos*(-RasActive + Total_Ras)/pow(KmSos - RasActive + Total_Ras, 2) - SosActive*cell*dkSos/(KmSos - RasActive + Total_Ras) + 2*SosActive*cell*kSos*(-RasActive + Total_Ras)*(-dKmSos + dRasActive - dTotal_Ras)/pow(KmSos - RasActive + Total_Ras, 3) + SosActive*cell*kSos*(-dRasActive + dTotal_Ras)/pow(KmSos - RasActive + Total_Ras, 2) - SosActive*cell*kSos*(-dKmSos + dRasActive - dTotal_Ras)/pow(KmSos - RasActive + Total_Ras, 2) + SosActive*dcell*kSos*(-RasActive + Total_Ras)/pow(KmSos - RasActive + Total_Ras, 2) - SosActive*dcell*kSos/(KmSos - RasActive + Total_Ras) - cell*dRasGapActive*kRasGap/(KmRasGap + RasActive) + cell*dSosActive*kSos*(-RasActive + Total_Ras)/pow(KmSos - RasActive + Total_Ras, 2) - cell*dSosActive*kSos/(KmSos - RasActive + Total_Ras);
pd[649] = cell*kSos*(-RasActive + Total_Ras)/(KmSos - RasActive + Total_Ras);
pd[651] = (RasActive*RasGapActive*cell*kRasGap/pow(KmRasGap + RasActive, 2) - RasGapActive*cell*kRasGap/(KmRasGap + RasActive) + SosActive*cell*kSos*(-RasActive + Total_Ras)/pow(KmSos - RasActive + Total_Ras, 2) - SosActive*cell*kSos/(KmSos - RasActive + Total_Ras)) + cj*(-1);
pd[666] = cell*dkRasToRaf1*(-Raf1Active + Total_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1) + cell*kRasToRaf1*(-Raf1Active + Total_Raf1)*(-dKmRasToRaf1 + dRaf1Active - dTotal_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 2) + cell*kRasToRaf1*(-dRaf1Active + dTotal_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1) + dcell*kRasToRaf1*(-Raf1Active + Total_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1);
pd[667] = AktActive*Raf1Active*cell*dkdRaf1ByAkt/pow(KmRaf1ByAkt + Raf1Active, 2) + 2*AktActive*Raf1Active*cell*kdRaf1ByAkt*(-dKmRaf1ByAkt - dRaf1Active)/pow(KmRaf1ByAkt + Raf1Active, 3) + AktActive*Raf1Active*dcell*kdRaf1ByAkt/pow(KmRaf1ByAkt + Raf1Active, 2) + AktActive*cell*dRaf1Active*kdRaf1ByAkt/pow(KmRaf1ByAkt + Raf1Active, 2) - AktActive*cell*dkdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) - AktActive*cell*kdRaf1ByAkt*(-dKmRaf1ByAkt - dRaf1Active)/pow(KmRaf1ByAkt + Raf1Active, 2) - AktActive*dcell*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) + Raf1Active*Raf1PPtase*cell*dkdRaf1/pow(KmdRaf1 + Raf1Active, 2) + 2*Raf1Active*Raf1PPtase*cell*kdRaf1*(-dKmdRaf1 - dRaf1Active)/pow(KmdRaf1 + Raf1Active, 3) + Raf1Active*Raf1PPtase*dcell*kdRaf1/pow(KmdRaf1 + Raf1Active, 2) + Raf1Active*cell*dAktActive*kdRaf1ByAkt/pow(KmRaf1ByAkt + Raf1Active, 2) + Raf1Active*cell*dRaf1PPtase*kdRaf1/pow(KmdRaf1 + Raf1Active, 2) + Raf1PPtase*cell*dRaf1Active*kdRaf1/pow(KmdRaf1 + Raf1Active, 2) - Raf1PPtase*cell*dkdRaf1/(KmdRaf1 + Raf1Active) - Raf1PPtase*cell*kdRaf1*(-dKmdRaf1 - dRaf1Active)/pow(KmdRaf1 + Raf1Active, 2) - Raf1PPtase*dcell*kdRaf1/(KmdRaf1 + Raf1Active) + RasActive*cell*dkRasToRaf1*(-Raf1Active + Total_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 2) - RasActive*cell*dkRasToRaf1/(KmRasToRaf1 - Raf1Active + Total_Raf1) + 2*RasActive*cell*kRasToRaf1*(-Raf1Active + Total_Raf1)*(-dKmRasToRaf1 + dRaf1Active - dTotal_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 3) + RasActive*cell*kRasToRaf1*(-dRaf1Active + dTotal_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 2) - RasActive*cell*kRasToRaf1*(-dKmRasToRaf1 + dRaf1Active - dTotal_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 2) + RasActive*dcell*kRasToRaf1*(-Raf1Active + Total_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 2) - RasActive*dcell*kRasToRaf1/(KmRasToRaf1 - Raf1Active + Total_Raf1) - cell*dAktActive*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) - cell*dRaf1PPtase*kdRaf1/(KmdRaf1 + Raf1Active) + cell*dRasActive*kRasToRaf1*(-Raf1Active + Total_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 2) - cell*dRasActive*kRasToRaf1/(KmRasToRaf1 - Raf1Active + Total_Raf1);
pd[672] = -Raf1Active*cell*dkdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) - Raf1Active*cell*kdRaf1ByAkt*(-dKmRaf1ByAkt - dRaf1Active)/pow(KmRaf1ByAkt + Raf1Active, 2) - Raf1Active*dcell*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) - cell*dRaf1Active*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active);
pd[681] = cell*kRasToRaf1*(-Raf1Active + Total_Raf1)/(KmRasToRaf1 - Raf1Active + Total_Raf1);
pd[682] = (AktActive*Raf1Active*cell*kdRaf1ByAkt/pow(KmRaf1ByAkt + Raf1Active, 2) - AktActive*cell*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active) + Raf1Active*Raf1PPtase*cell*kdRaf1/pow(KmdRaf1 + Raf1Active, 2) - Raf1PPtase*cell*kdRaf1/(KmdRaf1 + Raf1Active) + RasActive*cell*kRasToRaf1*(-Raf1Active + Total_Raf1)/pow(KmRasToRaf1 - Raf1Active + Total_Raf1, 2) - RasActive*cell*kRasToRaf1/(KmRasToRaf1 - Raf1Active + Total_Raf1)) + cj*(-1);
pd[687] = -Raf1Active*cell*kdRaf1ByAkt/(KmRaf1ByAkt + Raf1Active);
pd[698] = BRafActive*Raf1PPtase*cell*dkdBRaf/pow(BRafActive + KmdBRaf, 2) + 2*BRafActive*Raf1PPtase*cell*kdBRaf*(-dBRafActive - dKmdBRaf)/pow(BRafActive + KmdBRaf, 3) + BRafActive*Raf1PPtase*dcell*kdBRaf/pow(BRafActive + KmdBRaf, 2) + BRafActive*cell*dRaf1PPtase*kdBRaf/pow(BRafActive + KmdBRaf, 2) + Raf1PPtase*cell*dBRafActive*kdBRaf/pow(BRafActive + KmdBRaf, 2) - Raf1PPtase*cell*dkdBRaf/(BRafActive + KmdBRaf) - Raf1PPtase*cell*kdBRaf*(-dBRafActive - dKmdBRaf)/pow(BRafActive + KmdBRaf, 2) - Raf1PPtase*dcell*kdBRaf/(BRafActive + KmdBRaf) + Rap1Active*cell*dkRap1ToBRaf*(-BRafActive + Total_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 2) - Rap1Active*cell*dkRap1ToBRaf/(-BRafActive + KmRap1ToBRaf + Total_BRaf) + 2*Rap1Active*cell*kRap1ToBRaf*(-BRafActive + Total_BRaf)*(dBRafActive - dKmRap1ToBRaf - dTotal_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 3) + Rap1Active*cell*kRap1ToBRaf*(-dBRafActive + dTotal_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 2) - Rap1Active*cell*kRap1ToBRaf*(dBRafActive - dKmRap1ToBRaf - dTotal_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 2) + Rap1Active*dcell*kRap1ToBRaf*(-BRafActive + Total_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 2) - Rap1Active*dcell*kRap1ToBRaf/(-BRafActive + KmRap1ToBRaf + Total_BRaf) - cell*dRaf1PPtase*kdBRaf/(BRafActive + KmdBRaf) + cell*dRap1Active*kRap1ToBRaf*(-BRafActive + Total_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 2) - cell*dRap1Active*kRap1ToBRaf/(-BRafActive + KmRap1ToBRaf + Total_BRaf);
pd[704] = cell*dkRap1ToBRaf*(-BRafActive + Total_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf) + cell*kRap1ToBRaf*(-BRafActive + Total_BRaf)*(dBRafActive - dKmRap1ToBRaf - dTotal_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 2) + cell*kRap1ToBRaf*(-dBRafActive + dTotal_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf) + dcell*kRap1ToBRaf*(-BRafActive + Total_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf);
pd[713] = (BRafActive*Raf1PPtase*cell*kdBRaf/pow(BRafActive + KmdBRaf, 2) - Raf1PPtase*cell*kdBRaf/(BRafActive + KmdBRaf) + Rap1Active*cell*kRap1ToBRaf*(-BRafActive + Total_BRaf)/pow(-BRafActive + KmRap1ToBRaf + Total_BRaf, 2) - Rap1Active*cell*kRap1ToBRaf/(-BRafActive + KmRap1ToBRaf + Total_BRaf)) + cj*(-1);
pd[719] = cell*kRap1ToBRaf*(-BRafActive + Total_BRaf)/(-BRafActive + KmRap1ToBRaf + Total_BRaf);
pd[727] = cell*dkpRaf1*(-MekActive + Total_Mek)/(KmpRaf1 - MekActive + Total_Mek) + cell*kpRaf1*(-MekActive + Total_Mek)*(-dKmpRaf1 + dMekActive - dTotal_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 2) + cell*kpRaf1*(-dMekActive + dTotal_Mek)/(KmpRaf1 - MekActive + Total_Mek) + dcell*kpRaf1*(-MekActive + Total_Mek)/(KmpRaf1 - MekActive + Total_Mek);
pd[728] = cell*dkpBRaf*(-MekActive + Total_Mek)/(KmpBRaf - MekActive + Total_Mek) + cell*kpBRaf*(-MekActive + Total_Mek)*(-dKmpBRaf + dMekActive - dTotal_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 2) + cell*kpBRaf*(-dMekActive + dTotal_Mek)/(KmpBRaf - MekActive + Total_Mek) + dcell*kpBRaf*(-MekActive + Total_Mek)/(KmpBRaf - MekActive + Total_Mek);
pd[729] = BRafActive*cell*dkpBRaf*(-MekActive + Total_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 2) - BRafActive*cell*dkpBRaf/(KmpBRaf - MekActive + Total_Mek) + 2*BRafActive*cell*kpBRaf*(-MekActive + Total_Mek)*(-dKmpBRaf + dMekActive - dTotal_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 3) + BRafActive*cell*kpBRaf*(-dMekActive + dTotal_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 2) - BRafActive*cell*kpBRaf*(-dKmpBRaf + dMekActive - dTotal_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 2) + BRafActive*dcell*kpBRaf*(-MekActive + Total_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 2) - BRafActive*dcell*kpBRaf/(KmpBRaf - MekActive + Total_Mek) + MekActive*PP2AActive*cell*dkdMek/pow(KmdMek + MekActive, 2) + 2*MekActive*PP2AActive*cell*kdMek*(-dKmdMek - dMekActive)/pow(KmdMek + MekActive, 3) + MekActive*PP2AActive*dcell*kdMek/pow(KmdMek + MekActive, 2) + MekActive*cell*dPP2AActive*kdMek/pow(KmdMek + MekActive, 2) + PP2AActive*cell*dMekActive*kdMek/pow(KmdMek + MekActive, 2) - PP2AActive*cell*dkdMek/(KmdMek + MekActive) - PP2AActive*cell*kdMek*(-dKmdMek - dMekActive)/pow(KmdMek + MekActive, 2) - PP2AActive*dcell*kdMek/(KmdMek + MekActive) + Raf1Active*cell*dkpRaf1*(-MekActive + Total_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 2) - Raf1Active*cell*dkpRaf1/(KmpRaf1 - MekActive + Total_Mek) + 2*Raf1Active*cell*kpRaf1*(-MekActive + Total_Mek)*(-dKmpRaf1 + dMekActive - dTotal_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 3) + Raf1Active*cell*kpRaf1*(-dMekActive + dTotal_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 2) - Raf1Active*cell*kpRaf1*(-dKmpRaf1 + dMekActive - dTotal_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 2) + Raf1Active*dcell*kpRaf1*(-MekActive + Total_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 2) - Raf1Active*dcell*kpRaf1/(KmpRaf1 - MekActive + Total_Mek) + cell*dBRafActive*kpBRaf*(-MekActive + Total_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 2) - cell*dBRafActive*kpBRaf/(KmpBRaf - MekActive + Total_Mek) - cell*dPP2AActive*kdMek/(KmdMek + MekActive) + cell*dRaf1Active*kpRaf1*(-MekActive + Total_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 2) - cell*dRaf1Active*kpRaf1/(KmpRaf1 - MekActive + Total_Mek);
pd[742] = cell*kpRaf1*(-MekActive + Total_Mek)/(KmpRaf1 - MekActive + Total_Mek);
pd[743] = cell*kpBRaf*(-MekActive + Total_Mek)/(KmpBRaf - MekActive + Total_Mek);
pd[744] = (BRafActive*cell*kpBRaf*(-MekActive + Total_Mek)/pow(KmpBRaf - MekActive + Total_Mek, 2) - BRafActive*cell*kpBRaf/(KmpBRaf - MekActive + Total_Mek) + MekActive*PP2AActive*cell*kdMek/pow(KmdMek + MekActive, 2) - PP2AActive*cell*kdMek/(KmdMek + MekActive) + Raf1Active*cell*kpRaf1*(-MekActive + Total_Mek)/pow(KmpRaf1 - MekActive + Total_Mek, 2) - Raf1Active*cell*kpRaf1/(KmpRaf1 - MekActive + Total_Mek)) + cj*(-1);
pd[759] = cell*dkpMekCytoplasmic*(-ErkActive + Total_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk) + cell*kpMekCytoplasmic*(-ErkActive + Total_Erk)*(dErkActive - dKmpMekCytoplasmic - dTotal_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 2) + cell*kpMekCytoplasmic*(-dErkActive + dTotal_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk) + dcell*kpMekCytoplasmic*(-ErkActive + Total_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk);
pd[760] = ErkActive*PP2AActive*cell*dkdErk/pow(ErkActive + KmdErk, 2) + 2*ErkActive*PP2AActive*cell*kdErk*(-dErkActive - dKmdErk)/pow(ErkActive + KmdErk, 3) + ErkActive*PP2AActive*dcell*kdErk/pow(ErkActive + KmdErk, 2) + ErkActive*cell*dPP2AActive*kdErk/pow(ErkActive + KmdErk, 2) + MekActive*cell*dkpMekCytoplasmic*(-ErkActive + Total_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 2) - MekActive*cell*dkpMekCytoplasmic/(-ErkActive + KmpMekCytoplasmic + Total_Erk) + 2*MekActive*cell*kpMekCytoplasmic*(-ErkActive + Total_Erk)*(dErkActive - dKmpMekCytoplasmic - dTotal_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 3) + MekActive*cell*kpMekCytoplasmic*(-dErkActive + dTotal_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 2) - MekActive*cell*kpMekCytoplasmic*(dErkActive - dKmpMekCytoplasmic - dTotal_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 2) + MekActive*dcell*kpMekCytoplasmic*(-ErkActive + Total_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 2) - MekActive*dcell*kpMekCytoplasmic/(-ErkActive + KmpMekCytoplasmic + Total_Erk) + PP2AActive*cell*dErkActive*kdErk/pow(ErkActive + KmdErk, 2) - PP2AActive*cell*dkdErk/(ErkActive + KmdErk) - PP2AActive*cell*kdErk*(-dErkActive - dKmdErk)/pow(ErkActive + KmdErk, 2) - PP2AActive*dcell*kdErk/(ErkActive + KmdErk) + cell*dMekActive*kpMekCytoplasmic*(-ErkActive + Total_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 2) - cell*dMekActive*kpMekCytoplasmic/(-ErkActive + KmpMekCytoplasmic + Total_Erk) - cell*dPP2AActive*kdErk/(ErkActive + KmdErk);
pd[774] = cell*kpMekCytoplasmic*(-ErkActive + Total_Erk)/(-ErkActive + KmpMekCytoplasmic + Total_Erk);
pd[775] = (ErkActive*PP2AActive*cell*kdErk/pow(ErkActive + KmdErk, 2) + MekActive*cell*kpMekCytoplasmic*(-ErkActive + Total_Erk)/pow(-ErkActive + KmpMekCytoplasmic + Total_Erk, 2) - MekActive*cell*kpMekCytoplasmic/(-ErkActive + KmpMekCytoplasmic + Total_Erk) - PP2AActive*cell*kdErk/(ErkActive + KmdErk)) + cj*(-1);
pd[782] = cell*dkPI3K*(-PI3KActive + Total_PI3K)/(KmPI3K - PI3KActive + Total_PI3K) + cell*kPI3K*(-PI3KActive + Total_PI3K)*(-dKmPI3K + dPI3KActive - dTotal_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 2) + cell*kPI3K*(-dPI3KActive + dTotal_PI3K)/(KmPI3K - PI3KActive + Total_PI3K) + dcell*kPI3K*(-PI3KActive + Total_PI3K)/(KmPI3K - PI3KActive + Total_PI3K);
pd[786] = cell*dkPI3KRas*(-PI3KActive + Total_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K) + cell*kPI3KRas*(-PI3KActive + Total_PI3K)*(-dKmPI3KRas + dPI3KActive - dTotal_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 2) + cell*kPI3KRas*(-dPI3KActive + dTotal_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K) + dcell*kPI3KRas*(-PI3KActive + Total_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K);
pd[791] = RasActive*cell*dkPI3KRas*(-PI3KActive + Total_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 2) - RasActive*cell*dkPI3KRas/(KmPI3KRas - PI3KActive + Total_PI3K) + 2*RasActive*cell*kPI3KRas*(-PI3KActive + Total_PI3K)*(-dKmPI3KRas + dPI3KActive - dTotal_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 3) + RasActive*cell*kPI3KRas*(-dPI3KActive + dTotal_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 2) - RasActive*cell*kPI3KRas*(-dKmPI3KRas + dPI3KActive - dTotal_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 2) + RasActive*dcell*kPI3KRas*(-PI3KActive + Total_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 2) - RasActive*dcell*kPI3KRas/(KmPI3KRas - PI3KActive + Total_PI3K) + boundEGFReceptor*cell*dkPI3K*(-PI3KActive + Total_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 2) - boundEGFReceptor*cell*dkPI3K/(KmPI3K - PI3KActive + Total_PI3K) + 2*boundEGFReceptor*cell*kPI3K*(-PI3KActive + Total_PI3K)*(-dKmPI3K + dPI3KActive - dTotal_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 3) + boundEGFReceptor*cell*kPI3K*(-dPI3KActive + dTotal_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 2) - boundEGFReceptor*cell*kPI3K*(-dKmPI3K + dPI3KActive - dTotal_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 2) + boundEGFReceptor*dcell*kPI3K*(-PI3KActive + Total_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 2) - boundEGFReceptor*dcell*kPI3K/(KmPI3K - PI3KActive + Total_PI3K) + cell*dRasActive*kPI3KRas*(-PI3KActive + Total_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 2) - cell*dRasActive*kPI3KRas/(KmPI3KRas - PI3KActive + Total_PI3K) + cell*dboundEGFReceptor*kPI3K*(-PI3KActive + Total_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 2) - cell*dboundEGFReceptor*kPI3K/(KmPI3K - PI3KActive + Total_PI3K);
pd[797] = cell*kPI3K*(-PI3KActive + Total_PI3K)/(KmPI3K - PI3KActive + Total_PI3K);
pd[801] = cell*kPI3KRas*(-PI3KActive + Total_PI3K)/(KmPI3KRas - PI3KActive + Total_PI3K);
pd[806] = (RasActive*cell*kPI3KRas*(-PI3KActive + Total_PI3K)/pow(KmPI3KRas - PI3KActive + Total_PI3K, 2) - RasActive*cell*kPI3KRas/(KmPI3KRas - PI3KActive + Total_PI3K) + boundEGFReceptor*cell*kPI3K*(-PI3KActive + Total_PI3K)/pow(KmPI3K - PI3KActive + Total_PI3K, 2) - boundEGFReceptor*cell*kPI3K/(KmPI3K - PI3KActive + Total_PI3K)) + cj*(-1);
pd[821] = cell*dkAkt*(-AktActive + Total_Akt)/(-AktActive + KmAkt + Total_Akt) + cell*kAkt*(-AktActive + Total_Akt)*(dAktActive - dKmAkt - dTotal_Akt)/pow(-AktActive + KmAkt + Total_Akt, 2) + cell*kAkt*(-dAktActive + dTotal_Akt)/(-AktActive + KmAkt + Total_Akt) + dcell*kAkt*(-AktActive + Total_Akt)/(-AktActive + KmAkt + Total_Akt);
pd[822] = PI3KActive*cell*dkAkt*(-AktActive + Total_Akt)/pow(-AktActive + KmAkt + Total_Akt, 2) - PI3KActive*cell*dkAkt/(-AktActive + KmAkt + Total_Akt) + 2*PI3KActive*cell*kAkt*(-AktActive + Total_Akt)*(dAktActive - dKmAkt - dTotal_Akt)/pow(-AktActive + KmAkt + Total_Akt, 3) + PI3KActive*cell*kAkt*(-dAktActive + dTotal_Akt)/pow(-AktActive + KmAkt + Total_Akt, 2) - PI3KActive*cell*kAkt*(dAktActive - dKmAkt - dTotal_Akt)/pow(-AktActive + KmAkt + Total_Akt, 2) + PI3KActive*dcell*kAkt*(-AktActive + Total_Akt)/pow(-AktActive + KmAkt + Total_Akt, 2) - PI3KActive*dcell*kAkt/(-AktActive + KmAkt + Total_Akt) + cell*dPI3KActive*kAkt*(-AktActive + Total_Akt)/pow(-AktActive + KmAkt + Total_Akt, 2) - cell*dPI3KActive*kAkt/(-AktActive + KmAkt + Total_Akt);
pd[836] = cell*kAkt*(-AktActive + Total_Akt)/(-AktActive + KmAkt + Total_Akt);
pd[837] = (PI3KActive*cell*kAkt*(-AktActive + Total_Akt)/pow(-AktActive + KmAkt + Total_Akt, 2) - PI3KActive*cell*kAkt/(-AktActive + KmAkt + Total_Akt)) + cj*(-1);
pd[843] = cell*dkC3GNGF*(-C3GActive + Total_C3G)/(-C3GActive + KmC3GNGF + Total_C3G) + cell*kC3GNGF*(-C3GActive + Total_C3G)*(dC3GActive - dKmC3GNGF - dTotal_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 2) + cell*kC3GNGF*(-dC3GActive + dTotal_C3G)/(-C3GActive + KmC3GNGF + Total_C3G) + dcell*kC3GNGF*(-C3GActive + Total_C3G)/(-C3GActive + KmC3GNGF + Total_C3G);
pd[853] = boundNGFReceptor*cell*dkC3GNGF*(-C3GActive + Total_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 2) - boundNGFReceptor*cell*dkC3GNGF/(-C3GActive + KmC3GNGF + Total_C3G) + 2*boundNGFReceptor*cell*kC3GNGF*(-C3GActive + Total_C3G)*(dC3GActive - dKmC3GNGF - dTotal_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 3) + boundNGFReceptor*cell*kC3GNGF*(-dC3GActive + dTotal_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 2) - boundNGFReceptor*cell*kC3GNGF*(dC3GActive - dKmC3GNGF - dTotal_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 2) + boundNGFReceptor*dcell*kC3GNGF*(-C3GActive + Total_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 2) - boundNGFReceptor*dcell*kC3GNGF/(-C3GActive + KmC3GNGF + Total_C3G) + cell*dboundNGFReceptor*kC3GNGF*(-C3GActive + Total_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 2) - cell*dboundNGFReceptor*kC3GNGF/(-C3GActive + KmC3GNGF + Total_C3G);
pd[858] = cell*kC3GNGF*(-C3GActive + Total_C3G)/(-C3GActive + KmC3GNGF + Total_C3G);
pd[868] = (boundNGFReceptor*cell*kC3GNGF*(-C3GActive + Total_C3G)/pow(-C3GActive + KmC3GNGF + Total_C3G, 2) - boundNGFReceptor*cell*kC3GNGF/(-C3GActive + KmC3GNGF + Total_C3G)) + cj*(-1);
pd[883] = cell*dkC3G*(-Rap1Active + Total_Rap1)/(KmC3G - Rap1Active + Total_Rap1) + cell*kC3G*(-Rap1Active + Total_Rap1)*(-dKmC3G + dRap1Active - dTotal_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 2) + cell*kC3G*(-dRap1Active + dTotal_Rap1)/(KmC3G - Rap1Active + Total_Rap1) + dcell*kC3G*(-Rap1Active + Total_Rap1)/(KmC3G - Rap1Active + Total_Rap1);
pd[884] = C3GActive*cell*dkC3G*(-Rap1Active + Total_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 2) - C3GActive*cell*dkC3G/(KmC3G - Rap1Active + Total_Rap1) + 2*C3GActive*cell*kC3G*(-Rap1Active + Total_Rap1)*(-dKmC3G + dRap1Active - dTotal_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 3) + C3GActive*cell*kC3G*(-dRap1Active + dTotal_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 2) - C3GActive*cell*kC3G*(-dKmC3G + dRap1Active - dTotal_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 2) + C3GActive*dcell*kC3G*(-Rap1Active + Total_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 2) - C3GActive*dcell*kC3G/(KmC3G - Rap1Active + Total_Rap1) + Rap1Active*RapGapActive*cell*dkRapGap/pow(KmRapGap + Rap1Active, 2) + 2*Rap1Active*RapGapActive*cell*kRapGap*(-dKmRapGap - dRap1Active)/pow(KmRapGap + Rap1Active, 3) + Rap1Active*RapGapActive*dcell*kRapGap/pow(KmRapGap + Rap1Active, 2) + Rap1Active*cell*dRapGapActive*kRapGap/pow(KmRapGap + Rap1Active, 2) + RapGapActive*cell*dRap1Active*kRapGap/pow(KmRapGap + Rap1Active, 2) - RapGapActive*cell*dkRapGap/(KmRapGap + Rap1Active) - RapGapActive*cell*kRapGap*(-dKmRapGap - dRap1Active)/pow(KmRapGap + Rap1Active, 2) - RapGapActive*dcell*kRapGap/(KmRapGap + Rap1Active) + cell*dC3GActive*kC3G*(-Rap1Active + Total_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 2) - cell*dC3GActive*kC3G/(KmC3G - Rap1Active + Total_Rap1) - cell*dRapGapActive*kRapGap/(KmRapGap + Rap1Active);
pd[898] = cell*kC3G*(-Rap1Active + Total_Rap1)/(KmC3G - Rap1Active + Total_Rap1);
pd[899] = (C3GActive*cell*kC3G*(-Rap1Active + Total_Rap1)/pow(KmC3G - Rap1Active + Total_Rap1, 2) - C3GActive*cell*kC3G/(KmC3G - Rap1Active + Total_Rap1) + Rap1Active*RapGapActive*cell*kRapGap/pow(KmRapGap + Rap1Active, 2) - RapGapActive*cell*kRapGap/(KmRapGap + Rap1Active)) + cj*(-1);
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
double dcell = constants[69];
double dRasGapActive = constants[70];
double dRapGapActive = constants[71];
double dPP2AActive = constants[72];
double dRaf1PPtase = constants[73];
double dEGF_IC = constants[74];
double dNGF_IC = constants[75];
double dTotal_EGFReceptor = constants[76];
double dTotal_NGFReceptor = constants[77];
double dTotal_Sos = constants[78];
double dTotal_P90Rsk = constants[79];
double dTotal_Ras = constants[80];
double dTotal_Raf1 = constants[81];
double dTotal_BRaf = constants[82];
double dTotal_Mek = constants[83];
double dTotal_Erk = constants[84];
double dTotal_PI3K = constants[85];
double dTotal_Akt = constants[86];
double dTotal_C3G = constants[87];
double dTotal_Rap1 = constants[88];
double dZero = constants[89];
double dkrbEGF = constants[90];
double dkruEGF = constants[91];
double dkrbNGF = constants[92];
double dkruNGF = constants[93];
double dkEGF = constants[94];
double dKmEGF = constants[95];
double dkNGF = constants[96];
double dKmNGF = constants[97];
double dkdSos = constants[98];
double dKmdSos = constants[99];
double dkSos = constants[100];
double dKmSos = constants[101];
double dkRasGap = constants[102];
double dKmRasGap = constants[103];
double dkRasToRaf1 = constants[104];
double dKmRasToRaf1 = constants[105];
double dkpRaf1 = constants[106];
double dKmpRaf1 = constants[107];
double dkpBRaf = constants[108];
double dKmpBRaf = constants[109];
double dkdMek = constants[110];
double dKmdMek = constants[111];
double dkpMekCytoplasmic = constants[112];
double dKmpMekCytoplasmic = constants[113];
double dkdErk = constants[114];
double dKmdErk = constants[115];
double dkpP90Rsk = constants[116];
double dKmpP90Rsk = constants[117];
double dkPI3K = constants[118];
double dKmPI3K = constants[119];
double dkPI3KRas = constants[120];
double dKmPI3KRas = constants[121];
double dkAkt = constants[122];
double dKmAkt = constants[123];
double dkdRaf1ByAkt = constants[124];
double dKmRaf1ByAkt = constants[125];
double dkC3GNGF = constants[126];
double dKmC3GNGF = constants[127];
double dkC3G = constants[128];
double dKmC3G = constants[129];
double dkRapGap = constants[130];
double dKmRapGap = constants[131];
double dkRap1ToBRaf = constants[132];
double dKmRap1ToBRaf = constants[133];
double dkdRaf1 = constants[134];
double dKmdRaf1 = constants[135];
double dkdBRaf = constants[136];
double dKmdBRaf = constants[137];

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
dynamicVars[15] = dEGF_IC;
dynamicVars[16] = dNGF_IC;
dynamicVars[17] = dZero;
dynamicVars[18] = dZero;
dynamicVars[19] = dZero;
dynamicVars[20] = dZero;
dynamicVars[21] = dZero;
dynamicVars[22] = dZero;
dynamicVars[23] = dZero;
dynamicVars[24] = dZero;
dynamicVars[25] = dZero;
dynamicVars[26] = dZero;
dynamicVars[27] = dZero;
dynamicVars[28] = dZero;
dynamicVars[29] = dZero;
}
