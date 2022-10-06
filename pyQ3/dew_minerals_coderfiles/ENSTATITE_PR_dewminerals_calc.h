
const char *ENSTATITE_PR_dewminerals_identifier(void);
const char *ENSTATITE_PR_dewminerals_name(void);
const char *ENSTATITE_PR_dewminerals_formula(void);
const double ENSTATITE_PR_dewminerals_mw(void);
const double *ENSTATITE_PR_dewminerals_elements(void);

double ENSTATITE_PR_dewminerals_g(double T, double P);
double ENSTATITE_PR_dewminerals_dgdt(double T, double P);
double ENSTATITE_PR_dewminerals_dgdp(double T, double P);
double ENSTATITE_PR_dewminerals_d2gdt2(double T, double P);
double ENSTATITE_PR_dewminerals_d2gdtdp(double T, double P);
double ENSTATITE_PR_dewminerals_d2gdp2(double T, double P);
double ENSTATITE_PR_dewminerals_d3gdt3(double T, double P);
double ENSTATITE_PR_dewminerals_d3gdt2dp(double T, double P);
double ENSTATITE_PR_dewminerals_d3gdtdp2(double T, double P);
double ENSTATITE_PR_dewminerals_d3gdp3(double T, double P);

double ENSTATITE_PR_dewminerals_s(double T, double P);
double ENSTATITE_PR_dewminerals_v(double T, double P);
double ENSTATITE_PR_dewminerals_cv(double T, double P);
double ENSTATITE_PR_dewminerals_cp(double T, double P);
double ENSTATITE_PR_dewminerals_dcpdt(double T, double P);
double ENSTATITE_PR_dewminerals_alpha(double T, double P);
double ENSTATITE_PR_dewminerals_beta(double T, double P);
double ENSTATITE_PR_dewminerals_K(double T, double P);
double ENSTATITE_PR_dewminerals_Kp(double T, double P);

