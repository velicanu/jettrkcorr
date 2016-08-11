#import "draw3rdjet.C"

void all()
{
  for ( int i = 0 ; i < 2 ; i++) 
  {
    for(int j = -1 ; j < 5 ; ++j)
    {
      draw3rdjet(i,j,1);
    }
  }
}
