// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "MyNLP.hpp"
#include "IpJournalist.hpp"

#include <cassert>
#include <random>
#include <algorithm>    // std::max

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

using namespace Ipopt;

/* Constructor. */
MyNLP::MyNLP(Index dim, Index T, Number reg_param):
   dim(dim),
   T(T),
   reg_param(reg_param),
   min_x(1e-12),
   Xspace(NULL),
   betaspace(NULL),
   pispace(NULL),
   alphaspace(NULL),
   X(NULL),
   beta(NULL),
   alpha(NULL),
   pi(NULL),
   Xbeta(NULL),
   beta_eLR(NULL),
   jnlst(new Journalist())
{ 
   this->jnlst->AddFileJournal("Matrices", "matrices.out", J_INSUPPRESSIBLE);

   // Allocate memory for the parameters
   this->Xspace = new DenseGenMatrixSpace(dim, T);
   this->betaspace = new DenseGenMatrixSpace(dim, 1);
   this->pispace = new DenseGenMatrixSpace(T, 1);
   this->alphaspace = new DenseGenMatrixSpace(dim, dim);

   this->X        = Xspace->MakeNewDenseGenMatrix();
   this->beta     = betaspace->MakeNewDenseGenMatrix();
   this->alpha    = alphaspace->MakeNewDenseGenMatrix();
   this->pi       = pispace->MakeNewDenseGenMatrix();
   this->Xbeta    = Xspace->MakeNewDenseGenMatrix();
   this->beta_eLR = betaspace->MakeNewDenseGenMatrix();

   // Some sanity checks
   assert(this->dim == this->alphaspace->NRows());
   assert(this->alphaspace->NRows() == this->alphaspace->NCols());
   assert(this->dim == this->betaspace->NRows());
   assert(1 == this->betaspace->NCols());
   assert(this->dim == this->Xspace->NRows());
   assert(this->T == this->Xspace->NCols());

   // Initialize values of the parameters
   //InitializeParametersRandom();
   InitializeParametersSeries();

   #ifdef DEGUG
   this->X->Print(*jnlst, (EJournalLevel) J_INSUPPRESSIBLE, (EJournalCategory) J_DBG, "X", 0);
   this->beta->Print(*jnlst, (EJournalLevel) J_INSUPPRESSIBLE, (EJournalCategory) J_DBG, "beta", 0);
   this->pi->Print(*jnlst, (EJournalLevel) J_INSUPPRESSIBLE, (EJournalCategory) J_DBG, "pi", 0);
   #endif 

   // Build NLP problem and its data from the parameters
   AssembleDataNLP();

   #ifdef DEGUG
   this->Xbeta->Print(*jnlst, (EJournalLevel) J_INSUPPRESSIBLE, (EJournalCategory) J_DBG, "Xbeta", 0);
   this->alpha->Print(*jnlst, (EJournalLevel) J_INSUPPRESSIBLE, (EJournalCategory) J_DBG, "alpha", 0);
   this->beta_eLR->Print(*jnlst, (EJournalLevel) J_INSUPPRESSIBLE, (EJournalCategory) J_DBG, "beta_eLR", 0);
   #endif
}

void MyNLP::InitializeParametersRandom()
{
   std::default_random_engine generator;
   std::normal_distribution<double> distribution(0.0, 1.0);
   
   Index idx = 0;
   Number *X_values = this->X->Values(); //elements are stored one column after each other

   for (Index j = 0; j < this->T; j++)
   {
      for (Index i = 0; i < this->dim; i++)
      {
         X_values[idx++] = distribution(generator);
      }
   }

   Number *beta_values = this->beta->Values(); //elements are stored one column after each other
   for (Index i = 0; i < this->dim; i++)
   {
      beta_values[i] = distribution(generator);
   }

   Number *pi_values = this->pi->Values(); //elements are stored one column after each other
   for (Index i = 0; i < this->T; i++)
   {
      pi_values[i] = distribution(generator);
   }
}

void MyNLP::InitializeParametersSeries()
{
   
   Index idx = 0;
   Number *X_values = this->X->Values(); //elements are stored one column after each other

   for (Index j = 0; j < this->T; j++)
   {
      for (Index i = 0; i < this->dim; i++)
      {
         X_values[idx] = idx + 1.0;
         idx++;
      }
   }

   Number *beta_values = this->beta->Values(); //elements are stored one column after each other
   for (Index i = 0; i < this->dim; i++)
   {
      beta_values[i] = i + 1.0;
   }

   Number *pi_values = this->pi->Values(); //elements are stored one column after each other
   for (Index i = 0; i < this->T; i++)
   {
      pi_values[i] = i + 1.0;
   }
}

void MyNLP::AssembleDataNLP()
{
   
   // Compute low rank vector used in the outer product
   // Xbeta = X;
   // for i=1:T
   //    Xbeta(:,i) = X(:,i).*beta'*sqrt(1/T);
   // end
   Number *X_values = this->X->Values(); //elements are stored one column after each other
   Number *beta_values = this->beta->Values(); //elements are stored one column after each other
   Number *Xbeta_values = this->Xbeta->Values(); //elements are stored one column after each other
   for (Index j = 0; j < this->T; j++)
   {
      for (Index i = 0; i < this->dim; i++)
      {
         Xbeta_values[j*this->dim + i] = X_values[j*this->dim + i] * beta_values[i] * sqrt(1.0/this->T);
      }
   }

   // Compute the other constant vector in the objective
   // beta_eLR = (-2/T)*beta.*((pi.*X) * e)';
   Number *beta_eLR_values = this->beta_eLR->Values(); //elements are stored one column after each other
   Number *pi_values = this->pi->Values(); //elements are stored one column after each other

   for (Index i = 0; i < this->dim; i++)
   {
      beta_eLR_values[i] = 0;
   }
   for (Index j = 0; j < this->T; j++)
   {
      for (Index i = 0; i < this->dim; i++)
      {
         beta_eLR_values[i] += X_values[j*this->dim + i] * pi_values[j];
      }
      
   }
   for (Index i = 0; i < this->dim; i++)
   {
      beta_eLR_values[i] *= (-2.0/this->T) * beta_values[i];
   }


   // Compute the outer product
   // alpha = Xbeta*Xbeta';
   Number *alpha_values = this->alpha->Values(); //elements are stored one column after each other

   Index nrows = this->alphaspace->NRows();
   Index ncols = this->alphaspace->NCols();

   std::fill_n(alpha_values, nrows*ncols, 0.0);

   for (Index i = 0; i < this->dim; i++)
   {
      for (Index j = 0; j < this->dim; j++)
      {
         for (Index k = 0; k < this->T; k++)
         {
            //TODO: super inefficient, optimize the memory access pattern! Improve to iterate over continuous memory addresses.
            alpha_values[j*this->dim + i] += Xbeta_values[k*this->dim + i] * Xbeta_values[k*this->dim + j];
         }
      }
   }
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
   nnz_h_lag = n/2.0*(1+n); //consider only lower triangular part

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

   for (Index i = 0; i < n; i++)
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
   for (Index i = 0; i < n; i++)
   {
      x[i] = 0.1 + i/(double)n;
   }

   #ifdef DEGUG
   // Dump to the file
   DenseGenMatrix *X0     = this->betaspace->MakeNewDenseGenMatrix();
   Number *Xvals = X0->Values();
   for (Index i = 0; i < this->dim; i++)
   {
      Xvals[i] = x[i];
   }
   X0->Print(*(this->jnlst), (EJournalLevel) J_INSUPPRESSIBLE, (EJournalCategory) J_DBG, "X0", 0);
   #endif

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
   // fun=(beta_eLR+W*alpha)*W'+eps_C*sum(W.*(log(max(W,1e-12)))); where W is the unknown
   obj_value = 0.0;

   Number *alpha_values = this->alpha->Values(); //elements are stored one column after each other
   Number *beta_eLR_values = this->beta_eLR->Values(); //elements are stored one column after each other
   Number tmp;
   Index idx;

   for (Index i = 0; i < n; i++)
   {
      tmp = 0.0;
      for (Index j = 0; j < n; j++)
      {
         //tmp = alpha[i,j]*x[j]
         idx = this->dim*j + i;
         tmp += alpha_values[idx] * x[j];
      }
      // f = x[i]*(alpha[i,:]*x[:]) + beta_eLR[i]*x[i] + eps*x[i](log(max(x[i],thresh)))
      obj_value += x[i]*tmp + beta_eLR_values[i]*x[i] + this->reg_param*x[i]*log(std::max(x[i],this->min_x));
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
   // grad = beta+2*(alpha*W')'+eps_C.*(log(max(W,1e-12))+ones(1,d))

   Number *alpha_values = this->alpha->Values(); //elements are stored one column after each other
   Number *beta_eLR_values = this->beta_eLR->Values(); //elements are stored one column after each other

   Index idx;

   for (Index i = 0; i < n; i++)
   {
      grad_f[i] = 0.0;
      for (Index j = 0; j < n; j++)
      {
         grad_f[i] = 0.0;

         for (Index j = 0; j < n; j++)
         {
            idx = this->dim*j + i; //TODO: Improve to iterate over continuous memory addresses
            grad_f[i] += alpha_values[idx]*x[j];
         }

         grad_f[i] *= 2.0;
         grad_f[i] += beta_eLR_values[i] + this->reg_param*(1.0 + log(std::max(x[i],this->min_x)));
      }
   }

   #ifdef DEGUG
   // Dump to the file
   DenseGenMatrix *grad     = this->betaspace->MakeNewDenseGenMatrix();
   Number *gradvals = grad->Values();
   for (Index i = 0; i < this->dim; i++)
   {
      gradvals[i] = grad_f[i];
   }
   grad->Print(*(this->jnlst), (EJournalLevel) J_INSUPPRESSIBLE, (EJournalCategory) J_DBG, "grad", 0);
   #endif

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
   Number sum = 0.0;
   for (Index i = 0; i < n; i++)
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

      for (Index i = 0; i < n; i++)
      {
         // FORTRAN_STYLE indexing
         iRow[i] = 1;
         jCol[i] = i+1;
      }
   }
   else
   {
      // return the values of the jacobian of the constraints
      for (Index i = 0; i < n; i++)
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
      // Ipopt expects that only the lower diagonal entries are specified

      Index nnz = 0;
      for (Index i = 0; i < n; i++)
      {
         for (Index j = 0; j <= i; j++)
         {
            // FORTRAN_STYLE indexing
            iRow[nnz] = i+1;
            jCol[nnz] = j+1;
            nnz++;
         }
      }
      assert(nnz == nele_hess);
   }
   else
   {
      // return the values of the Hessian of Lagrangian
      Number *alpha_values = this->alpha->Values(); //elements are stored one column after each other
      Index alpha_idx;
      Index nnz_idx = 0;
      
      for (Index i = 0; i < n; i++)
      {
         for (Index j = 0; j <= i; j++)
         {
            alpha_idx = this->dim*j + i;
            values[nnz_idx] = 2.0 * obj_factor * alpha_values[alpha_idx];

            // Add the regularization eps_C./max(x(i),1e-12)
            if (i == j)
            {
               values[nnz_idx] += obj_factor*this->reg_param/std::max(x[i], this->min_x);
            }

            //printf("ALPHA: %d %d %f\n",i,j, alpha_values[alpha_idx] );
            //printf("HESS(%f):  %d %d %f\n",obj_factor, i,j, values[nnz_idx] );

            nnz_idx++;
         }
      }
      assert(nnz_idx == nele_hess);

      #ifdef DEGUG
      // Dump to the file
      DenseGenMatrixSpace hessspace(nele_hess, 1);
      DenseGenMatrix *hess     = hessspace.MakeNewDenseGenMatrix();
      Number *hessvals = hess->Values();
      for (Index i = 0; i < nele_hess; i++)
      {
         hessvals[i] = values[i];
      }
      hess->Print(*(this->jnlst), (EJournalLevel) J_INSUPPRESSIBLE, (EJournalCategory) J_DBG, "hessian", 0);
      #endif
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
