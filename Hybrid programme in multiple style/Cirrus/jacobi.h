void jacobistep(double **psinew, double **psi, int m, int n, int myid);

void jacobistepvort(double **zetnew, double **psinew,
		    double **zet,    double** psi,
		    int m, int n, double re, int myid,int ln);

double deltasq(double **newarr, double **oldarr, int m, int n);
