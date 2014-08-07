{
double getG (int a, int b)
{
  double sumG;
  sumG = 0;
  for (int gi=0; gi<p.size(); gi++)
  {
    sumG += A[a][gi]*A[gi][b];
  }
  return sumG;
}

double getD (int a, int b)
{
  double sumD;
  sumD = 0;
  if (a==b)
  {
    for (int di=0; di<p.size(); di++)
    {
      return 1/sqrt(getG(a,b));
    }    
  }
  else
  {
    return 0;
  }
}

double getR (int a, int b)
{
  double sumR;
  sumR = 0;
  for (int t=0; t<p.size(); t++)
  {
    for (int k=0; k<p.size(); k++)
    {
      sumR+=(getD(a,k)*getG(k,t))*getD(t,b)

    }
  }
  return sumR;
}

}