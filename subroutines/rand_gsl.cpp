#include "gsl/gsl_rng.h"
#include "rand_gsl.h"

static gsl_rng *the_generator = NULL;

void srand_gsl48( int num)
{

  if( !the_generator )
    the_generator = gsl_rng_alloc( gsl_rng_ranlxd1 );
  gsl_rng_set( the_generator, num);

}

double rand_gsl48( void )
{
  if( !the_generator )
    the_generator = gsl_rng_alloc(gsl_rng_ranlxd1 );
  // if(!the_generator )
//     {
//       printf("Failed to load gsl random number generator.\n");
//       exit(1);
//     }
  return gsl_rng_uniform( the_generator );
}

void srand_gsl( int num)
{

  if( !the_generator )
    the_generator = gsl_rng_alloc( gsl_rng_taus );
  gsl_rng_set( the_generator, num);

}

double rand_gsl( void )
{
  if( !the_generator )
    the_generator = gsl_rng_alloc(gsl_rng_taus );
  // if(!the_generator )
//     {
//       printf("Failed to load gsl random number generator.\n");
//       exit(1);
//     }
  return gsl_rng_uniform( the_generator );
}
