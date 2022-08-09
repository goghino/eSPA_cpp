// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "MyNLP.hpp"

#include <cassert>
#include <random>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

using namespace Ipopt;

/* Constructor. */
MyNLP::MyNLP(int dim, int T, double reg_param):
   dim(dim),
   T(T),
   reg_param(reg_param),
   Xspace(NULL),
   betaspace(NULL),
   X(NULL),
   beta(NULL),
   pi(NULL),
   Xbeta(NULL),
   beta_eLR(NULL)
{ 
   // Allocate memory for the parameters
   this->Xspace = new DenseGenMatrixSpace(dim, T);
   this->betaspace = new DenseGenMatrixSpace(dim, 1);

   this->X        = Xspace->MakeNewDenseGenMatrix();
   this->beta     = betaspace->MakeNewDenseGenMatrix();
   this->pi       = betaspace->MakeNewDenseGenMatrix();
   this->Xbeta    = Xspace->MakeNewDenseGenMatrix();
   this->beta_eLR = betaspace->MakeNewDenseGenMatrix();

   // Initialize values of the parameters
   InitializeParameters();

   // Build NLP problem and its data from the parameters
   AssembleDataNLP();
}

void MyNLP::InitializeParameters()
{
   std::default_random_engine generator;
   std::normal_distribution<double> distribution(0.0, 1.0);
   
   int idx = 0;
   Number *X_values = this->X->Values(); //elements are stored one column after each other

   for (int j = 0; j < this->T; j++)
   {
      for (int i = 0; i < this->dim; i++)
      {
         X_values[idx++] = distribution(generator);
      }
   }


   Number *beta_values = this->beta->Values(); //elements are stored one column after each other
   Number *pi_values = this->pi->Values(); //elements are stored one column after each other
   for (int i = 0; i < this->dim; i++)
   {
      beta_values[i] = distribution(generator);
      pi_values[i] = distribution(generator);
   }
}

void MyNLP::AssembleDataNLP()
{
   
   // // Compute low rank vector used in the outer product
   // // Xbeta = X;
   // // for i=1:T
   // //    Xbeta(:,i) = X(:,i).*beta'*sqrt(1/T);
   // // end
   // for (int i = 0; i < this->dim; i++)
   // {
   //    for (int j = 0; j < this->T; j++)
   //    {
   //       this->Xbeta[i*this->T + j] = this->X[i*this->T + j] * beta[j] * sqrt(1/T);
   //    }
   // }

   // // Compute the other constant vector in the objective
   // // beta_eLR = (-2/T)*beta.*((pi.*X) * e)';
   // for (int i = 0; i < this->dim; i++)
   // {
   //    this->beta_eLR[i] = 0;
   //    for (int j = 0; j < this->T; j++)
   //    {
   //       this->beta_eLR[i] += this->X[i*this->T + j] * this->pi[j];
   //    }
   //    this->beta_eLR[i] *= (-2/this->T) * this->beta[i];
   // }
}

MyNLP::~MyNLP()
{ 
}

bool MyNLP::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
   // The problem described in MyNLP.hpp has 2 variables, x1, & x2,
   n = this->dim;

   // one equality constraint,
   m = 1;

   // nonzeros in the jacobian,
   nnz_jac_g = n;

   // nonzeros in the hessian of the lagrangian
   nnz_h_lag = (n+m)*(n+m);

   // We use the standard fortran index style for row/col entries
   index_style = FORTRAN_STYLE;

   return true;
}

bool MyNLP::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
   Number* g_l,
   Number* g_u
)
{
   // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
   // If desired, we could assert to make sure they are what we think they are.

   for (int i = 0; i < n; i++)
   {
      x_l[i] = 0.0;
      x_u[i] = 1.0;
   }

   // we have one equality constraint, so we set the bounds on this constraint
   // to be equal (and zero).
   g_l[0] = g_u[0] = 0.0;

   return true;
}

bool MyNLP::get_starting_point(
   Index   n,
   bool    init_x,
   Number* x,
   bool    init_z,
   Number* z_L,
   Number* z_U,
   Index   m,
   bool    init_lambda,
   Number* lambda
)
{
   // Here, we assume we only have starting values for x, if you code
   // your own NLP, you can provide starting values for the others if
   // you wish.
   assert(init_x == true);
   assert(init_z == false);
   assert(init_lambda == false);

   std::default_random_engine generator;
   std::normal_distribution<double> distribution(0.0, 1.0);

   // we initialize x in bounds, in the upper right quadrant
   for (int i = 0; i < n; i++)
   {
      x[i] = distribution(generator);
      printf("%d: %f\n",i, x[i]);
   }

   return true;
}

bool MyNLP::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
   // return the value of the objective function
   obj_value = 0.0;
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
      {
         obj_value += x[i]*x[j];
      }
   }

   return true;
}

bool MyNLP::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
   // return the gradient of the objective function grad_{x} f(x)
   for (int i = 0; i < n; i++)
   {
      grad_f[i] = 0.0;
      for (int j = 0; j < n; j++)
      {
         grad_f[i] += x[j];
      }
   }

   return true;
}

bool MyNLP::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
   // return the value of the constraints: g(x)
   double sum = 0.0;
   for (int i = 0; i < n; i++)
   {
      sum += x[i];
   }

   g[0] = sum - 1;

   return true;
}

bool MyNLP::eval_jac_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Index         nele_jac,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
   if( values == NULL )
   {
      // return the structure of the jacobian of the constraints

      for (int i = 0; i < n; i++)
      {
         // FORTRAN_STYLE indexing
         iRow[i] = 1;
         jCol[i] = i+1;
      }
   }
   else
   {
      // return the values of the jacobian of the constraints
      for (int i = 0; i < n; i++)
      {
         // i-th element
         values[i] = 1.0;
      }
   }

   return true;
}

bool MyNLP::eval_h(
   Index         n,
   const Number* x,
   bool          new_x,
   Number        obj_factor,
   Index         m,
   const Number* lambda,
   bool          new_lambda,
   Index         nele_hess,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
   if( values == NULL )
   {
      // return the structure of the Hessian of Lagrangian

      int nnz = 0;
      for (int i = 0; i < n+m; i++)
      {
         for (int j = 0; j < n+m; j++)
         {
            if (i == n+m && j == n+m)
            {
               continue;
            }

            // FORTRAN_STYLE indexing
            iRow[nnz] = i+1;
            jCol[nnz] = j+1;
            nnz++;
         }
      }
   }
   else
   {
      // return the values of the Hessian of Lagrangian
      for (int i = 0; i < nele_hess; i++)
      {
         // i-th element
         values[i] = 1.0 * obj_factor;
      }
   }

   return true;
}

void MyNLP::finalize_solution(
   SolverReturn               status,
   Index                      n,
   const Number*              x,
   const Number*              z_L,
   const Number*              z_U,
   Index                      m,
   const Number*              g,
   const Number*              lambda,
   Number                     obj_value,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
)
{
   // here is where we would store the solution to variables, or write to a file, etc
   // so we could use the solution. Since the solution is displayed to the console,
   // we currently do nothing here.
}
