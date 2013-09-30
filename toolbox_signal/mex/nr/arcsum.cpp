void arcsum(unsigned long iin[], unsigned long iout[], unsigned long ja,
	int nwk, unsigned long nrad, unsigned long nc)
{
	int j,karry=0;
	unsigned long jtmp;

	for (j=nwk;j>nc;j--) {
		jtmp=ja;
		ja /= nrad;
		iout[j]=iin[j]+(jtmp-ja*nrad)+karry;
		if (iout[j] >= nrad) {
			iout[j] -= nrad;
			karry=1;
		} else karry=0;
	}
	iout[nc]=iin[nc]+ja+karry;
}
