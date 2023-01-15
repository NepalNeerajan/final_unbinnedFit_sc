Double_t funcA(double *x, double *p){
    double a = 0;
    a = p[0]*exp(-p[1]*x[0]);
    return a*p[1];
}

Double_t funcB(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a4 = exp(-p[4]*x[0]);
    a = p[0]*(1-p[2]-p[3])*p[1]*(a1/(p[4]-p[1]) + a4/(p[1]-p[4]));
    return a*p[4];
}

Double_t funcC(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a4 = exp(-p[4]*x[0]);
    double a7 = exp(-p[7]*x[0]);
    a = p[0]*(1-p[2]-p[3])*(1-p[5]-p[6])*p[1]*p[4]*(a1/(p[4]-p[1])/(p[7]-p[1]) + a4/(p[1]-p[4])/(p[7]-p[4]) + a7/(p[1]-p[7])/(p[4]-p[7]));
    return a*p[7];
}

Double_t funcD(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a4 = exp(-p[4]*x[0]);
    double a7 = exp(-p[7]*x[0]);
    double a10 = exp(-p[10]*x[0]);
    a = p[0]*(1-p[2]-p[3])*(1-p[5]-p[6])*(1-p[8]-p[9])*p[1]*p[4]*p[7]*(a1/(p[4]-p[1])/(p[7]-p[1])/(p[10]-p[1]) + a4/(p[1]-p[4])/(p[7]-p[4])/(p[10]-p[4]) + a7/(p[1]-p[7])/(p[4]-p[7])/(p[10]-p[7]) + a10/(p[1]-p[10])/(p[7]-p[10])/(p[4]-p[10]));
    return a*p[10];
}

Double_t funcE(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a13 = exp(-p[13]*x[0]);
    a = p[0]*p[2]*p[1]*(a1/(p[13]-p[1]) + a13/(p[1]-p[13]));
    return a*p[13];
}

Double_t funcF(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a4 = exp(-p[4]*x[0]);
    double a13 = exp(-p[13]*x[0]);
    double a16 = exp(-p[16]*x[0]);
    double A_B_F = p[0]*(1-p[2]-p[3])*p[5]*p[1]*p[4]*(a1/(p[4]-p[1])/(p[16]-p[1]) + a4/(p[1]-p[4])/(p[16]-p[4]) + a16/(p[1]-p[16])/(p[4]-p[16]));
    double A_E_F = p[0]*p[2]*(1-p[14]-p[15])*p[1]*p[13]*(a1/(p[13]-p[1])/(p[16]-p[1]) + a13/(p[1]-p[13])/(p[16]-p[13]) + a16/(p[1]-p[16])/(p[13]-p[16]));
    a = A_B_F + A_E_F;
    return a*p[16];
}

Double_t funcG(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a4 = exp(-p[4]*x[0]);
    double a7 = exp(-p[7]*x[0]);
    double a13 = exp(-p[13]*x[0]);
    double a16 = exp(-p[16]*x[0]);
    double a19 = exp(-p[19]*x[0]);
    double A_B_C_G = p[0]*(1-p[2]-p[3])*(1-p[5]-p[6])*p[8]*p[1]*p[4]*p[7]*(a1/(p[4]-p[1])/(p[7]-p[1])/(p[19]-p[1]) + a4/(p[1]-p[4])/(p[7]-p[4])/(p[19]-p[4]) + a7/(p[4]-p[7])/(p[1]-p[7])/(p[19]-p[7]) + a19/(p[1]-p[19])/(p[4]-p[19])/(p[7]-p[19]));
    double A_B_F_G = p[0]*(1-p[2]-p[3])*p[5]*(1-p[17]-p[18])*p[1]*p[4]*p[16]*(a1/(p[4]-p[1])/(p[16]-p[1])/(p[19]-p[1]) + a4/(p[1]-p[4])/(p[16]-p[4])/(p[19]-p[4]) + a16/(p[4]-p[16])/(p[1]-p[16])/(p[19]-p[16]) + a19/(p[1]-p[19])/(p[4]-p[19])/(p[16]-p[19]));
    double A_E_F_G = p[0]*p[2]*(1-p[14]-p[15])*(1-p[17]-p[18])*p[1]*p[13]*p[16]*(a1/(p[13]-p[1])/(p[16]-p[1])/(p[19]-p[1]) + a13/(p[1]-p[13])/(p[16]-p[13])/(p[19]-p[13]) + a16/(p[13]-p[16])/(p[1]-p[16])/(p[19]-p[16]) + a19/(p[1]-p[19])/(p[13]-p[19])/(p[16]-p[19]));
    a = A_B_C_G + A_B_F_G + A_E_F_G;
    return a*p[19];
}


Double_t funcH(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a22 = exp(-p[22]*x[0]);
    a = p[0]*p[3]*p[1]*(a1/(p[22]-p[1]) + a22/(p[1]-p[22]));
    return a*p[22];
}

Double_t funcI(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a4 = exp(-p[4]*x[0]);
    double a13 = exp(-p[13]*x[0]);
    double a22 = exp(-p[22]*x[0]);
    double a25 = exp(-p[25]*x[0]);
    double A_B_I = p[0]*(1-p[2]-p[3])*p[6]*p[1]*p[4]*(a1/(p[4]-p[1])/(p[25]-p[1]) + a4/(p[1]-p[4])/(p[25]-p[4]) + a25/(p[1]-p[25])/(p[4]-p[25]));
    double A_E_I = p[0]*p[2]*p[14]*p[1]*p[13]*(a1/(p[13]-p[1])/(p[25]-p[1]) + a13/(p[1]-p[13])/(p[25]-p[13]) + a25/(p[1]-p[25])/(p[13]-p[25]));
    double A_H_I = p[0]*p[3]*(1-p[23]-p[24])*p[1]*p[22]*(a1/(p[22]-p[1])/(p[25]-p[1]) + a22/(p[1]-p[22])/(p[25]-p[22]) + a25/(p[1]-p[25])/(p[22]-p[25]));
    a = A_B_I + A_E_I + A_H_I;
    return a*p[25];
}

Double_t funcJ(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a4 = exp(-p[4]*x[0]);
    double a7 = exp(-p[7]*x[0]);
    double a13 = exp(-p[13]*x[0]);
    double a16 = exp(-p[16]*x[0]);
    double a22 = exp(-p[22]*x[0]);
    double a25 = exp(-p[25]*x[0]);
    double a28 = exp(-p[28]*x[0]);
    double A_B_C_J = p[0]*(1-p[2]-p[3])*(1-p[5]-p[6])*p[9]*p[1]*p[4]*p[7]*(a1/(p[4]-p[1])/(p[7]-p[1])/(p[28]-p[1]) + a4/(p[1]-p[4])/(p[7]-p[4])/(p[28]-p[4]) + a7/(p[4]-p[7])/(p[1]-p[7])/(p[28]-p[7]) + a28/(p[4]-p[28])/(p[7]-p[28])/(p[1]-p[28]));
    double A_B_F_J = p[0]*(1-p[2]-p[3])*p[5]*p[17]*p[1]*p[4]*p[16]*(a1/(p[4]-p[1])/(p[16]-p[1])/(p[28]-p[1]) + a4/(p[1]-p[4])/(p[16]-p[4])/(p[28]-p[4]) + a16/(p[4]-p[16])/(p[1]-p[16])/(p[28]-p[16]) + a28/(p[4]-p[28])/(p[16]-p[28])/(p[1]-p[28]));
    double A_B_I_J = p[0]*(1-p[2]-p[3])*p[6]*(1-p[26]-p[27])*p[1]*p[4]*p[25]*(a1/(p[4]-p[1])/(p[25]-p[1])/(p[28]-p[1]) + a4/(p[1]-p[4])/(p[25]-p[4])/(p[28]-p[4]) + a25/(p[4]-p[25])/(p[1]-p[25])/(p[28]-p[25]) + a28/(p[4]-p[28])/(p[25]-p[28])/(p[1]-p[28]));
    double A_E_F_J = p[0]*p[2]*(1-p[14]-p[15])*p[17]*p[1]*p[13]*p[16]*(a1/(p[13]-p[1])/(p[16]-p[1])/(p[28]-p[1]) + a13/(p[1]-p[13])/(p[16]-p[13])/(p[28]-p[13]) + a16/(p[13]-p[16])/(p[1]-p[16])/(p[28]-p[16]) + a28/(p[13]-p[28])/(p[16]-p[28])/(p[1]-p[28]));
    double A_E_I_J = p[0]*p[2]*p[14]*(1-p[26]-p[27])*p[1]*p[13]*p[25]*(a1/(p[13]-p[1])/(p[25]-p[1])/(p[28]-p[1]) + a13/(p[1]-p[13])/(p[25]-p[13])/(p[28]-p[13]) + a25/(p[13]-p[25])/(p[1]-p[25])/(p[28]-p[25]) + a28/(p[13]-p[28])/(p[25]-p[28])/(p[1]-p[28]));
    double A_H_I_J = p[0]*p[3]*(1-p[23]-p[24])*(1-p[26]-p[27])*p[1]*p[22]*p[25]*(a1/(p[22]-p[1])/(p[25]-p[1])/(p[28]-p[1]) + a22/(p[1]-p[22])/(p[25]-p[22])/(p[28]-p[22]) + a25/(p[22]-p[25])/(p[1]-p[25])/(p[28]-p[25]) + a28/(p[22]-p[28])/(p[25]-p[28])/(p[1]-p[28]));
    a = A_B_C_J + A_B_F_J + A_B_I_J + A_E_F_J + A_E_I_J + A_H_I_J;
    return a*p[28];
}

Double_t funcK(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a13 = exp(-p[13]*x[0]);
    double a22 = exp(-p[22]*x[0]);
    double a31 = exp(-p[31]*x[0]);
    double A_E_K = p[0]*p[2]*p[15]*p[1]*p[13]*(a1/(p[13]-p[1])/(p[31]-p[1]) + a13/(p[1]-p[13])/(p[31]-p[13]) + a31/(p[13]-p[31])/(p[1]-p[31]));
    double A_H_K = p[0]*p[3]*p[23]*p[1]*p[22]*(a1/(p[22]-p[1])/(p[31]-p[1]) + a22/(p[1]-p[22])/(p[31]-p[22]) + a31/(p[22]-p[31])/(p[1]-p[31]));
    a = A_E_K + A_H_K;
    return a*p[31];
}

Double_t funcL(double *x, double *p){
    double a = 0;
    double a1 = exp(-p[1]*x[0]);
    double a4 = exp(-p[4]*x[0]);
    double a13 = exp(-p[13]*x[0]);
    double a16 = exp(-p[16]*x[0]);
    double a22 = exp(-p[22]*x[0]);
    double a25 = exp(-p[25]*x[0]);
    double a31 = exp(-p[31]*x[0]);
    double a34 = exp(-p[34]*x[0]);
    double A_B_F_L = p[0]*(1-p[2]-p[3])*p[5]*p[18]*p[1]*p[4]*p[16]*(a1/(p[4]-p[1])/(p[16]-p[1])/(p[34]-p[1]) + a4/(p[1]-p[4])/(p[16]-p[4])/(p[34]-p[4]) + a16/(p[4]-p[16])/(p[1]-p[16])/(p[34]-p[16]) + a34/(p[1]-p[34])/(p[16]-p[34])/(p[4]-p[34]));
    double A_B_I_L = p[0]*(1-p[2]-p[3])*p[6]*p[26]*p[1]*p[4]*p[25]*(a1/(p[4]-p[1])/(p[25]-p[1])/(p[34]-p[1]) + a4/(p[1]-p[4])/(p[25]-p[4])/(p[34]-p[4]) + a25/(p[4]-p[25])/(p[1]-p[25])/(p[34]-p[25]) + a34/(p[1]-p[34])/(p[25]-p[34])/(p[4]-p[34]));
    double A_E_F_L = p[0]*p[2]*(1-p[14]-p[15])*p[18]*p[1]*p[13]*p[16]*(a1/(p[13]-p[1])/(p[16]-p[1])/(p[34]-p[1]) + a13/(p[1]-p[13])/(p[16]-p[13])/(p[34]-p[13]) + a16/(p[13]-p[16])/(p[1]-p[16])/(p[34]-p[16]) + a34/(p[1]-p[34])/(p[16]-p[34])/(p[13]-p[34]));
    double A_E_I_L = p[0]*p[2]*p[14]*p[26]*p[1]*p[13]*p[25]*(a1/(p[13]-p[1])/(p[25]-p[1])/(p[34]-p[1]) + a13/(p[1]-p[13])/(p[25]-p[13])/(p[34]-p[13]) + a25/(p[13]-p[25])/(p[1]-p[25])/(p[34]-p[25]) + a34/(p[1]-p[34])/(p[25]-p[34])/(p[13]-p[34]));
    double A_H_I_L = p[0]*p[3]*(1-p[23]-p[24])*p[26]*p[1]*p[22]*p[25]*(a1/(p[22]-p[1])/(p[25]-p[1])/(p[34]-p[1]) + a22/(p[1]-p[22])/(p[25]-p[22])/(p[34]-p[22]) + a25/(p[22]-p[25])/(p[1]-p[25])/(p[34]-p[25]) + a34/(p[1]-p[34])/(p[25]-p[34])/(p[22]-p[34]));
    double A_H_K_L = p[0]*p[3]*p[23]*(1-p[32]-p[33])*p[1]*p[22]*p[31]*(a1/(p[22]-p[1])/(p[31]-p[1])/(p[34]-p[1]) + a22/(p[1]-p[22])/(p[31]-p[22])/(p[34]-p[22]) + a31/(p[22]-p[31])/(p[1]-p[31])/(p[34]-p[31]) + a34/(p[1]-p[34])/(p[31]-p[34])/(p[22]-p[34]));
    a = A_B_F_L + A_B_I_L + A_E_F_L + A_E_I_L + A_H_I_L + A_H_K_L;
    return a*p[34];
}

//to draw the individual isotopes with added bkg
Double_t funcA_d(double *x, double *p){
	return funcA(x,p)+p[36];
}

Double_t funcB_d(double *x, double *p){
	return funcB(x,p) + p[36];
}

Double_t funcC_d(double *x, double *p){
	return funcC(x,p) + p[36];
}


Double_t funcE_d(double *x, double *p){
	return funcE(x,p) + p[36];
}

Double_t funcH_d(double *x, double *p){
	return funcH(x,p) + p[36];
}

Double_t func36_d(double *x, double *p){
	return p[36];
}






//defining the global function
Double_t FT(double *x, double *p){
    double returnVal = p[36];
    if(x[0]>-5 && x[0]<5){
       TF1::RejectPoint();
       return returnVal;
     }
    double func = funcA(x,p) + funcB(x,p) + funcC(x,p) + funcD(x,p) + funcE(x,p) + funcF(x,p) + funcG(x,p) + funcH(x,p) + funcI(x,p) + funcJ(x,p) + funcK(x,p) + funcL(x,p);
    if(x[0]>0.0) return returnVal += func;
    return returnVal;
}

Double_t F1NT(double *x, double *p){
    double returnVal = p[37];
    if(x[0]>-5 && x[0]<5){
       TF1::RejectPoint();
       return returnVal;
     }
    double func = (funcA(x,p)*p[2] + funcB(x,p)*p[5] + funcC(x,p)*p[8] + funcD(x,p)*p[11] + funcE(x,p)*p[14] + funcF(x,p)*p[17] + funcG(x,p)*p[20] + funcH(x,p)*p[23] + funcI(x,p)*p[26] + funcJ(x,p)*p[29])*p[35] + 2*p[41]*(1-p[41])*(funcA(x,p)*p[3] + funcB(x,p)*p[6] + funcC(x,p)*p[9] + funcD(x,p)*p[12] + funcE(x,p)*p[15] + funcF(x,p)*p[18] + funcG(x,p)*p[21] + funcH(x,p)*p[24]) + p[39]*FT(x,p);
    if(x[0]>0.0) return returnVal += func;
    return returnVal;

}

Double_t F2NT(double *x, double *p){
    double returnVal = p[38];
    if(x[0]>-5 && x[0]<5){
       TF1::RejectPoint();
       return returnVal;
     }
    double func = p[41]*p[41]*(funcA(x,p)*p[3] + funcB(x,p)*p[6] + funcC(x,p)*p[9] + funcD(x,p)*p[12] + funcE(x,p)*p[15] + funcF(x,p)*p[18] + funcG(x,p)*p[21] +funcH(x,p)*p[24]) + p[40]*FT(x,p) + p[39]*F1NT(x,p);
    if(x[0]>0.0) return returnVal += func;
    return returnVal;
}



//to draw in the P1n plot
Double_t funcA_d1(double *x, double *p){ //d
  return funcA(x,p)*p[2]*p[35] + 2*p[41]*(1-p[41])*funcA(x,p)*p[3] + p[37];
}

Double_t funcB_d1(double *x, double *p){ //d
	return funcB(x,p)*p[5]*p[35] + 2*p[41]*(1-p[41])*funcB(x,p)*p[6] + p[37];
}

Double_t funcC_d1(double *x, double *p){ //d
	return funcC(x,p)*p[8]*p[35] +  2*p[41]*(1-p[41])*funcC(x,p)*p[9] + p[37];
}


Double_t funcE_d1(double *x, double *p){ //d
  return funcE(x,p)*p[14]*p[35] +  2*p[41]*(1-p[41])*funcE(x,p)*p[15] + p[37];
}


Double_t funcH_d1(double *x, double *p){ //d
  return funcD(x,p)*p[23]*p[35] + 2*p[41]*(1-p[41])*funcH(x,p)*p[24] + p[37];
}



//for decay 2n
Double_t funcA_d2(double *x, double *p){ //d
  return funcA(x,p)*p[3]*p[41]*p[41] + p[38];
}

Double_t funcB_d2(double *x, double *p){ //d
	return funcB(x,p)*p[6]*p[41]*p[41]  +  p[38];
}

Double_t funcC_d2(double *x, double *p){ //d
	return funcC(x,p)*p[9]*p[41]*p[41]  +  p[38];
}


Double_t funcE_d2(double *x, double *p){ //d
	return funcE(x,p)*p[15]*p[41]*p[41] +  p[38];
}

Double_t funcH_d2(double *x, double *p){ //d
	return funcH(x,p)*p[24]*p[41]*p[41] +  p[38];
}

