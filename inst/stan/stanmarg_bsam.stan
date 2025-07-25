/* This file is based on LERSIL.stan by Ben Goodrich.
   https://github.com/bgoodri/LERSIL */
functions { // you can use these in R following `rstan::expose_stan_functions("foo.stan")`
  /*
    Fills in the elements of a coefficient matrix containing some mix of 
    totally free, free subject to a sign constraint, and fixed elements
    
    @param free_elements vector of unconstrained elements
    @param skeleton matrix of the same dimensions as the output whose elements are
      positive_infinity(): if output element is totally free
      other: if output element is fixed to that number
    @return matrix of coefficients
  */
  matrix fill_matrix(vector free_elements, matrix skeleton, array[,] int eq_skeleton, int pos_start, int spos_start) {
    int R = rows(skeleton);
    int C = cols(skeleton);
    matrix[R, C] out;

    int pos = spos_start; // position of eq_skeleton
    int freepos = pos_start; // position of free_elements
    int eqelem = 0;
    
    for (c in 1:C) for (r in 1:R) {
      real rc = skeleton[r, c];
      if (is_inf(rc)) { // free
	int eq = eq_skeleton[pos, 1];
	int wig = eq_skeleton[pos, 3];
	if (eq == 0 || wig == 1) {
	  out[r,c] = free_elements[freepos];
	  freepos += 1;
	} else {
	  eqelem = eq_skeleton[pos, 2];
	  out[r,c] = free_elements[eqelem];
	}
	pos += 1;
      } else out[r,c] = skeleton[r, c]; // fixed, so do not bump pos
    }
    return out;
  }

  /*
    Constructs vector of unconstrained elements/parameters, based on filled
    parameter matrix. (this is the reverse of fill_matrix())
 */
  vector fill_vector(array[] matrix parm_matrix, array[] matrix skeleton, array[,] int eq_skeleton, int len_free, int Ng) {
    int R = rows(skeleton[1]);
    int C = cols(skeleton[1]);
    vector[len_free] out;

    int pos = 1; // position of eq_skeleton
    int freepos = 1; // position of free_elements

    for (g in 1:Ng) {
      for (c in 1:C) for (r in 1:R) {
        real rc = skeleton[g, r, c];
        if (is_inf(rc)) { // free
          int eq = eq_skeleton[pos, 1];
	  int wig = eq_skeleton[pos, 3];
	  if (eq == 0 || wig == 1) {
	    out[freepos] = parm_matrix[g,r,c];
	    freepos += 1;
	  }
	  pos += 1;
        }
      }
    }
    return out;
  }
  
  vector fill_prior(vector free_elements, array[] real pri_mean, array[,] int eq_skeleton) {
    int R = dims(eq_skeleton)[1];
    int eqelem = 0;
    int pos = 1;
    vector[num_elements(pri_mean)] out;

    for (r in 1:R) {
      if (pos <= num_elements(pri_mean)) {
	int eq = eq_skeleton[r, 1];
	int wig = eq_skeleton[r, 3];

	if (eq == 0) {
	  out[pos] = pri_mean[pos];
	  pos += 1;
	} else if (wig == 1) {
	  eqelem = eq_skeleton[r, 2];
	  out[pos] = free_elements[eqelem];
	  pos += 1;
	}
      }
    }
    return out;
  }
  
  /*
   * This is a bug-free version of csr_to_dense_matrix and has the same arguments
   */
  matrix to_dense_matrix(int m, int n, vector w, array[] int v, array[] int u) {
    matrix[m, n] out = rep_matrix(0, m, n);
    int pos = 1;
    for (i in 1:m) {
      int start = u[i];
      int nnz = u[i + 1] - start;
      for (j in 1:nnz) {
        out[i, v[pos]] = w[pos];
        pos += 1;
      }
    }
    return out;
  }

  // sign function
  int sign(real x) {
    if (x > 0)
      return 1;
    else
      return -1;
  }

  // sign-constrain a vector of loadings
  vector sign_constrain_load(vector free_elements, int npar, array[,] int sign_mat) {
    vector[npar] out;
    for (i in 1:npar) {
      if (sign_mat[i,1]) {
        int lookupval = sign_mat[i,2];
        if (free_elements[lookupval] < 0) {
	  out[i] = -free_elements[i];
	} else {
	  out[i] = free_elements[i];
	}
      } else {
        out[i] = free_elements[i];
      }
    }
    return out;
  }

  // sign-constrain a vector of regressions or covariances
  vector sign_constrain_reg(vector free_elements, int npar, array[,] int sign_mat, vector load_par1, vector load_par2) {
    vector[npar] out;
    for (i in 1:npar) {
      if (sign_mat[i,1]) {
        int lookupval1 = sign_mat[i,2];
	int lookupval2 = sign_mat[i,3];
        if (sign(load_par1[lookupval1]) * sign(load_par2[lookupval2]) < 0) {
	  out[i] = -free_elements[i];
	} else {
	  out[i] = free_elements[i];
	}
      } else {
        out[i] = free_elements[i];
      }
    }
    return out;
  }

  // obtain covariance parameter vector for correlation/sd matrices
  vector cor2cov(array[] matrix cormat, array[] matrix sdmat, int num_free_elements, array[] matrix matskel, array[,] int wskel, int ngrp) {
    vector[num_free_elements] out;
    int R = rows(to_matrix(cormat[1]));
    int pos = 1; // position of eq_skeleton
    int freepos = 1; // position of free_elements
    
    for (g in 1:ngrp) {
      for (c in 1:(R-1)) for (r in (c+1):R) {
        if (is_inf(matskel[g,r,c])) {
	  if (wskel[pos,1] == 0) {
	    out[freepos] = sdmat[g,r,r] * sdmat[g,c,c] * cormat[g,r,c];
	    freepos += 1;
	  }
	  pos += 1;
	}
      }
    }
    return out;
  }

  // E step of EM algorithm on latent continuous space
  array[] matrix estep(array[] vector YXstar, array[] vector Mu, array[] matrix Sigma, array[] int Nobs, array[,] int Obsvar, array[] int startrow, array[] int endrow, array[] int grpnum, int Np, int Ng) {
    int p = dims(YXstar)[2];
    array[Ng] matrix[p, p + 1] out; //mean vec + cov mat
    matrix[dims(YXstar)[1], p] YXfull; // columns consistenly ordered
    matrix[p, p] T2pat;
    array[p] int obsidx;
    int r1;
    int r2;
    int grpidx;
    int Nmis;

    for (g in 1:Ng) {
      out[g] = rep_matrix(0, p, p + 1);
    }

    for (mm in 1:Np) {
      obsidx = Obsvar[mm,];
      r1 = startrow[mm];
      r2 = endrow[mm];
      grpidx = grpnum[mm];
      Nmis = p - Nobs[mm];

      if (Nobs[mm] < p) {
	matrix[Nobs[mm], Nobs[mm]] Sig22 = Sigma[grpidx, obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]];
	matrix[Nmis, Nmis] Sig11 = Sigma[grpidx, obsidx[(Nobs[mm] + 1):p], obsidx[(Nobs[mm] + 1):p]];
	matrix[Nmis, Nobs[mm]] Sig12 = Sigma[grpidx, obsidx[(Nobs[mm] + 1):p], obsidx[1:Nobs[mm]]];
	matrix[Nobs[mm], Nobs[mm]] S22inv = inverse_spd(Sig22);
	matrix[Nmis, Nmis] T2p11 = Sig11 - (Sig12 * S22inv * Sig12');
	
        // partition into observed/missing, compute Sigmas, add to out
	for (jj in r1:r2) {
	  vector[Nmis] ymis;
	  ymis = Mu[grpidx, obsidx[(Nobs[mm] + 1):p]] + (Sig12 * S22inv * (YXstar[jj, 1:Nobs[mm]] - Mu[grpidx, obsidx[1:Nobs[mm]]]));
	  for (kk in 1:Nobs[mm]) {
	    YXfull[jj, obsidx[kk]] = YXstar[jj, kk];
	  }
	  for (kk in (Nobs[mm] + 1):p) {
	    YXfull[jj, obsidx[kk]] = ymis[kk - Nobs[mm]];
	  }
	}
	T2pat = crossprod(YXfull[r1:r2,]);
	// correction for missing cells/conditional covariances
	for (jj in 1:Nmis) {
	  for (kk in jj:Nmis) {
	    T2pat[obsidx[Nobs[mm] + jj], obsidx[Nobs[mm] + kk]] = T2pat[obsidx[Nobs[mm] + jj], obsidx[Nobs[mm] + kk]] + (r2 - r1 + 1) * T2p11[jj, kk];
	    if (kk > jj) {
	      T2pat[obsidx[Nobs[mm] + kk], obsidx[Nobs[mm] + jj]] = T2pat[obsidx[Nobs[mm] + jj], obsidx[Nobs[mm] + kk]];
	    }
	  }
	}
      } else {
	// complete data
	for (jj in r1:r2) {
	  for (kk in 1:Nobs[mm]) {
	    YXfull[jj, obsidx[kk]] = YXstar[jj, kk];
	  }
	}
	T2pat = crossprod(YXfull[r1:r2,]);
      }
      for (i in 1:p) {
	out[grpidx,i,1] += sum(YXfull[r1:r2,i]);
      }
      out[grpidx,,2:(p+1)] += T2pat;
    }
    
    return out;
  }

  matrix sig_inv_update(matrix Sigmainv, array[] int obsidx, int Nobs, int np, real logdet) {
    matrix[Nobs + 1, Nobs + 1] out = rep_matrix(0, Nobs + 1, Nobs + 1);
    int nrm = np - Nobs;
    matrix[nrm, nrm] H;
    matrix[nrm, Nobs] A;

    if (nrm == 0) {
      out[1:Nobs, 1:Nobs] = Sigmainv;
      out[Nobs + 1, Nobs + 1] = logdet;
    } else {
      H = Sigmainv[obsidx[(Nobs + 1):np], obsidx[(Nobs + 1):np]];
      A = Sigmainv[obsidx[(Nobs + 1):np], obsidx[1:Nobs]];

      out[1:Nobs, 1:Nobs] = Sigmainv[obsidx[1:Nobs], obsidx[1:Nobs]] - A' * mdivide_left_spd(H, A);
      out[Nobs + 1, Nobs + 1] = logdet + log_determinant(H);
    }

    return out;
  }
  
  real multi_normal_suff(vector xbar, matrix S, vector Mu, matrix Supdate, int N) {
    int Nobs = dims(S)[1];
    real out;

    // using elementwise multiplication + sum here for efficiency
    out = -.5 * N * ( sum(Supdate[1:Nobs, 1:Nobs] .* (S + (xbar - Mu) * (xbar - Mu)')) + Supdate[Nobs + 1, Nobs + 1] + Nobs * log(2 * pi()) );

    if(is_nan(out) || out == positive_infinity()) out = negative_infinity();

    return out;
  }

  real cond_density(real theta, vector Mu, vector SDvec, vector Lambda_y, matrix Tau, int nvar) {
    vector[nvar] condmn;
    array[nvar] real probs;
    real out;

    condmn = Mu + Lambda_y * theta;
    
    for (i in 1:nvar) {
      //array[2] real phis;

      //phis[1] = normal_lcdf(runtau[1] | condmn[i], SDvec[i]);
      //phis[2] = normal_lcdf(runtau[2] | condmn[i], SDvec[i]);

      //logprobs[i] = phis[2] + log1m_exp(phis[1] - phis[2]);
      probs[i] = bernoulli_lpmf(1 | Phi_approx((Tau[i,2] - condmn[i])/SDvec[i]) - Phi_approx((Tau[i,1] - condmn[i])/SDvec[i]));
      //if (is_inf(probs[i])) probs[i] = -20;
    }

    out = sum(probs);
    
    return out;
  }

  real trunc_normal_rng(real Mu, real SD, row_vector runtau) {
    real p_lb = normal_cdf(runtau[1] | Mu, SD);
    real p_ub = normal_cdf(runtau[2] | Mu, SD);
    real u = uniform_rng(p_lb, p_ub);
    real y = Mu + SD * inv_Phi(u);
    return y;
  }

  // fill covariance matrix with blocks
  array[] matrix fill_cov(array[] matrix covmat, array[,] int blkse, array[] int nblk,
			  array[] matrix mat_1, array[] matrix mat_2, array[] matrix mat_3,
			  array[] matrix mat_4, array[] matrix mat_5, array[,] int psiorder,
			  array[,] int psirevord) {
    int Ng = dims(covmat)[1];
    array[Ng] matrix[dims(covmat)[2], dims(covmat)[3]] out;

    for (g in 1:Ng) {
      out[g] = covmat[g, psiorder[g], psiorder[g]];
    }

    for (k in 1:sum(nblk)) {
      int blkidx = blkse[k, 6];
      int arrayidx = blkse[k, 5];
      int blkgrp = blkse[k, 4];
      int srow = blkse[k, 1];
      int erow = blkse[k, 2];

      if (arrayidx == 1) {
	out[blkgrp, srow:erow, srow:erow] = mat_1[blkidx];
      } else if (arrayidx == 2) {
	out[blkgrp, srow:erow, srow:erow] = mat_2[blkidx];
      } else if (arrayidx == 3) {
	out[blkgrp, srow:erow, srow:erow] = mat_3[blkidx];
      } else if (arrayidx == 4) {
	out[blkgrp, srow:erow, srow:erow] = mat_4[blkidx];
      } else {
	out[blkgrp, srow:erow, srow:erow] = mat_5[blkidx];
      }
    }

    for (g in 1:Ng) {
      out[g] = out[g, psirevord[g], psirevord[g]];
    }

    return out;
  }
}
data {
  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA2#v=onepage&q=LISREL&f=false
  int<lower=0> p; // number of manifest response variables
  int<lower=0> p_c; // number of manifest level 2 variables
  int<lower=0> q; // number of manifest predictors
  int<lower=0> m; // number of latent endogenous variables
  int<lower=0> m_c; // number of latent level 2 variables
  int<lower=0> n; // number of latent exogenous variables
  int<lower=1> Ng; // number of groups
  int<lower=0, upper=1> missing; // are there missing values?
  int<lower=0, upper=1> save_lvs; // should we save lvs?
  int<lower=1> Np; // number of group-by-missing patterns combos
  array[Ng] int<lower=1> N; // number of observations per group
  array[Np] int<lower=1> Nobs; // number of observed variables in each missing pattern
  array[Np] int<lower=0> Nordobs; // number of ordinal observed variables in each missing pattern
  array[Np, p + q] int<lower=0> Obsvar; // indexing of observed variables
  int<lower=1> Ntot; // number of observations across all groups
  array[Np] int<lower=1> startrow; // starting row for each missing pattern
  array[Np] int<lower=1,upper=Ntot> endrow; // ending row for each missing pattern
  array[Np] int<lower=1,upper=Ng> grpnum; // group number for each row of data
  int<lower=0,upper=1> wigind; // do any parameters have approx equality constraint ('wiggle')?
  int<lower=0, upper=1> has_data; // are the raw data on y and x available?
  int<lower=0, upper=1> ord; // are there any ordinal variables?
  int<lower=0, upper=1> multilev; // is this a multilevel dataset?
  int<lower=0> Nord; // how many ordinal variables?
  array[Nord] int<lower=0> ordidx; // indexing of ordinal variables
  array[Np, Nord] int<lower=0> OrdObsvar; // indexing of observed ordinal variables in YXo
  int<lower=0> Noent; // how many observed entries of ordinal variables (for data augmentation)
  array[p + q - Nord] int<lower=0> contidx; // indexing of continuous variables
  array[Nord] int<lower=1> nlevs; // how many levels does each ordinal variable have
  array[Ng, 2] int<lower=1> nclus; // number of level 1 + level 2 observations
  int<lower=0> p_tilde; // total number of variables
  array[Ntot] vector[multilev ? p_tilde : p + q - Nord] YX; // continuous data
  array[Ntot, Nord] int YXo; // ordinal data
  array[Np] int<lower=0> Nx; // number of fixed.x variables (within)
  array[Np] int<lower=0> Nx_between; // number of fixed.x variables (between)
  int<lower=0, upper=1> use_cov;
  int<lower=0, upper=1> pri_only;
  int<lower=0> emiter; // number of em iterations for saturated model in ppp (missing data only)
  int<lower=0, upper=1> use_suff; // should we compute likelihood via mvn sufficient stats?
  int<lower=0, upper=1> do_test; // should we do everything in generated quantities?
  array[Np] vector[multilev ? p_tilde : p + q - Nord] YXbar; // sample means of continuous manifest variables
  array[Np] matrix[multilev ? (p_tilde + 1) : (p + q - Nord + 1), multilev ? (p_tilde + 1) : (p + q - Nord + 1)] S;     // sample covariance matrix among all continuous manifest variables NB!! multiply by (N-1) to use wishart lpdf!!
  int<lower=1> ngibbs; // number of gibbs iterations for structural parameters per iteration
  array[Np] int<lower=0> Ndum; // number of ovs with dummy lvs
  array[Np, p] int<lower=1> dum_ov_idx; // first Ndum are dummy ovs, then non-dummy ovs
  array[Np, m] int<lower=1> dum_lv_idx; // first Ndum are dummy lvs corresponding to ovs, then non-dummy lvs
  array[Np] int<lower=0> Ndum_x; // number of eXo dummy ovs/lvx
  array[Np, p] int<lower=1> dum_ov_x_idx; // index of eXo dummy ovs/lvs
  array[Np, m] int<lower=1> dum_lv_x_idx;

  array[Ng] int<lower=1> measnblk; // number of blocks in measurement model
  array[sum(measnblk), 3] int<lower=0> measblkse; // start/end rows of blocks
  array[Ng, p] int<lower=1> measorder; // reordering to get blocks
  array[Ng, p] int<lower=1> measrevord; // reverse ordering
  
  int ngh;
  vector[ngh] ghnode;
  vector[ngh] ghwt;  
  
  array[sum(nclus[,2])] int<lower=1> cluster_size; // number of obs per cluster
  array[Ng] int<lower=1> ncluster_sizes; // number of unique cluster sizes
  array[sum(ncluster_sizes)] int<lower=1> cluster_sizes; // unique cluster sizes
  array[sum(ncluster_sizes)] int<lower=1> cluster_size_ns; // number of clusters of each size
  array[Np, multilev ? p_tilde : p + q] int<lower=0> Xvar; // indexing of fixed.x variables (within)
  array[Np, multilev ? p_tilde : p + q] int<lower=0> Xdatvar; // indexing of fixed.x in data (differs from Xvar when missing)
  array[Np, multilev ? p_tilde : p + q] int<lower=0> Xbetvar; // indexing of fixed.x variables (between)
  array[sum(ncluster_sizes)] vector[p_tilde] mean_d; // sample means by unique cluster size
  array[sum(ncluster_sizes)] matrix[p_tilde, p_tilde] cov_d; // sample covariances by unique cluster size
  array[Ng] matrix[p_tilde, p_tilde] cov_w; // observed "within" covariance matrix
  array[sum(nclus[,2])] vector[p_tilde] mean_d_full; // sample means/covs by cluster, for clusterwise log-densities
  array[sum(nclus[,2])] matrix[p_tilde, p_tilde] cov_d_full;
  array[Ng] vector[p_tilde] xbar_w; // data estimates of within/between means/covs (for saturated logl)
  array[Ng] vector[p_tilde] xbar_b;
  array[Ng] matrix[p_tilde, p_tilde] cov_b;
  array[Ng] real gs; // group size constant, for computation of saturated logl
  int N_within; // number of within variables
  int N_between; // number of between variables
  int N_both; // number of variables at both levels
  array[2] int N_lev; // number of observed variables at each level
  array[N_within] int within_idx;
  array[p_tilde] int between_idx; // between indexing, followed by within/both
  array[N_lev[1]] int ov_idx1;
  array[N_lev[2]] int ov_idx2;
  array[N_both] int both_idx;
  vector[multilev ? sum(ncluster_sizes) : Ng] log_lik_x; // ll of fixed x variables by unique cluster size
  vector[multilev ? sum(nclus[,2]) : Ng] log_lik_x_full; // ll of fixed x variables by cluster
  
  
  /* sparse matrix representations of skeletons of coefficient matrices, 
     which is not that interesting but necessary because you cannot pass
     missing values into the data block of a Stan program from R */
  int<lower=0> len_w1;        // max number of free elements in Lambda_y per grp
  array[Ng] int<lower=0> wg1;           // number of free elements in Lambda_y per grp
  array[Ng] vector[len_w1] w1;          // values of free elements in Lambda_y
  array[Ng, len_w1] int<lower=1> v1;    // index  of free elements in Lambda_y
  array[Ng, p + 1] int<lower=1> u1;     // index  of free elements in Lambda_y
  array[sum(wg1), 3] int<lower=0> w1skel;
  array[sum(wg1), 2] int<lower=0> lam_y_sign;
  int<lower=0> len_lam_y;     // number of free elements minus equality constraints
  array[len_lam_y] real lambda_y_mn;           // prior
  array[len_lam_y] real<lower=0> lambda_y_sd;

  // same things but for B
  int<lower=0> len_w4;
  array[Ng] int<lower=0> wg4;
  array[Ng] vector[len_w4] w4;
  array[Ng, len_w4] int<lower=1> v4;
  array[Ng, m + 1] int<lower=1> u4;
  array[sum(wg4), 3] int<lower=0> w4skel;
  array[sum(wg4), 3] int<lower=0> b_sign;
  int<lower=0> len_b;
  array[len_b] real b_mn;
  array[len_b] real<lower=0> b_sd;
  
  // same things but for diag(Theta)
  int<lower=0> len_w5;
  array[Ng] int<lower=0> wg5;
  array[Ng] vector[len_w5] w5;
  array[Ng, len_w5] int<lower=1> v5;
  array[Ng, p + 1] int<lower=1> u5;
  array[sum(wg5), 3] int<lower=0> w5skel;
  int<lower=0> len_thet_sd;
  array[len_thet_sd] real<lower=0> theta_sd_shape;
  array[len_thet_sd] real<lower=0> theta_sd_rate;
  int<lower=-2, upper=2> theta_pow;

  // same things but for Theta_r
  int<lower=0> len_w7;
  array[Ng] int<lower=0> wg7;
  array[Ng] vector[len_w7] w7;
  array[Ng, len_w7] int<lower=1> v7;
  array[Ng, p + 1] int<lower=1> u7;
  array[sum(wg7), 3] int<lower=0> w7skel;
  int<lower=0> len_thet_r;
  array[len_thet_r] real<lower=0> theta_r_alpha;
  array[len_thet_r] real<lower=0> theta_r_beta;

  // for blocks within Theta_r that receive lkj
  array[5] int<lower=0> thetanblk;
  array[5] int<lower=2> thetadims;
  array[sum(thetanblk), 7] int<lower=0> thetablkse;
  int<lower=0> len_w8;
  array[Ng] int<lower=0> wg8;
  array[Ng] vector[len_w8] w8;
  array[Ng, len_w8] int<lower=1> v8;
  array[Ng, p + 1] int<lower=1> u8;
  array[sum(wg8), 3] int<lower=0> w8skel;
  array[Ng, p] int<lower=1> thetaorder;
  array[Ng, p] int<lower=1> thetarevord;
  
  // same things but for Psi
  int<lower=0> len_w9;
  array[Ng] int<lower=0> wg9;
  array[Ng] vector[len_w9] w9;
  array[Ng, len_w9] int<lower=1> v9;
  array[Ng, m + 1] int<lower=1> u9;
  array[sum(wg9), 3] int<lower=0> w9skel;
  int<lower=0> len_psi_sd;
  array[len_psi_sd] real<lower=0> psi_sd_shape;
  array[len_psi_sd] real<lower=0> psi_sd_rate;
  int<lower=-2,upper=2> psi_pow;

  // same things but for Psi_r
  int<lower=0> len_w10;
  array[Ng] int<lower=0> wg10;
  array[Ng] vector[len_w10] w10;
  array[Ng, len_w10] int<lower=1> v10;
  array[Ng, m + 1] int<lower=1> u10;
  array[sum(wg10), 3] int<lower=0> w10skel;
  array[sum(wg10), 3] int<lower=0> psi_r_sign;
  int<lower=0> len_psi_r;
  array[len_psi_r] real<lower=0> psi_r_alpha;
  array[len_psi_r] real<lower=0> psi_r_beta;

  // for blocks within Psi_r
  array[5] int<lower=0> psinblk;
  array[5] int<lower=1> psidims;
  array[sum(psinblk), 7] int<lower=0> psiblkse;
  int<lower=0> len_w11;
  array[Ng] int<lower=0> wg11;
  array[Ng] vector[len_w11] w11;
  array[Ng, len_w11] int<lower=1> v11;
  array[Ng, m + 1] int<lower=1> u11;
  array[sum(wg11), 3] int<lower=0> w11skel;
  array[Ng, m] int<lower=1> psiorder;
  array[Ng, m] int<lower=1> psirevord;
  
  // same things but for Nu
  int<lower=0> len_w13;
  array[Ng] int<lower=0> wg13;
  array[Ng] vector[len_w13] w13;
  array[Ng, len_w13] int<lower=1> v13;
  array[Ng, use_cov ? 1 : p + q + 1] int<lower=1> u13;
  array[sum(wg13), 3] int<lower=0> w13skel;
  int<lower=0> len_nu;
  array[len_nu] real nu_mn;
  array[len_nu] real<lower=0> nu_sd;
  
  // same things but for Alpha
  int<lower=0> len_w14;
  array[Ng] int<lower=0> wg14;
  array[Ng] vector[len_w14] w14;
  array[Ng, len_w14] int<lower=0> v14;
  array[Ng, use_cov ? 1 : m + n + 1] int<lower=1> u14;
  array[sum(wg14), 3] int<lower=0> w14skel;
  array[sum(wg14), 3] int<lower=0> alph_sign;
  int<lower=0> len_alph;
  array[len_alph] real alpha_mn;
  array[len_alph] real<lower=0> alpha_sd;

  // same things but for Tau
  int<lower=0> len_w15;
  array[Ng] int<lower=0> wg15;
  array[Ng] vector[len_w15] w15;
  array[Ng, len_w15] int<lower=0> v15;
  array[Ng, sum(nlevs) - Nord + 1] int<lower=1> u15;
  array[sum(wg15), 3] int<lower=0> w15skel;
  int<lower=0> len_tau;
  array[len_tau] real tau_mn;
  array[len_tau] real<lower=0> tau_sd;

  // Level 2 matrices start here!!
  // Lambda
  int<lower=0> len_w1_c;
  array[Ng] int<lower=0> wg1_c;
  array[Ng] vector[len_w1_c] w1_c;
  array[Ng, len_w1_c] int<lower=1> v1_c;
  array[Ng, p_c + 1] int<lower=1> u1_c;
  array[sum(wg1_c), 3] int<lower=0> w1skel_c;
  array[sum(wg1_c), 2] int<lower=0> lam_y_sign_c;
  int<lower=0> len_lam_y_c;
  array[len_lam_y_c] real lambda_y_mn_c;
  array[len_lam_y_c] real<lower=0> lambda_y_sd_c;

  // same things but for B
  int<lower=0> len_w4_c;
  array[Ng] int<lower=0> wg4_c;
  array[Ng] vector[len_w4_c] w4_c;
  array[Ng, len_w4_c] int<lower=1> v4_c;
  array[Ng, m_c + 1] int<lower=1> u4_c;
  array[sum(wg4_c), 3] int<lower=0> w4skel_c;
  array[sum(wg4_c), 3] int<lower=0> b_sign_c;
  int<lower=0> len_b_c;
  array[len_b_c] real b_mn_c;
  array[len_b_c] real<lower=0> b_sd_c;
  
  // same things but for diag(Theta)
  int<lower=0> len_w5_c;
  array[Ng] int<lower=0> wg5_c;
  array[Ng] vector[len_w5_c] w5_c;
  array[Ng, len_w5_c] int<lower=1> v5_c;
  array[Ng, p_c + 1] int<lower=1> u5_c;
  array[sum(wg5_c), 3] int<lower=0> w5skel_c;
  int<lower=0> len_thet_sd_c;
  array[len_thet_sd_c] real<lower=0> theta_sd_shape_c;
  array[len_thet_sd_c] real<lower=0> theta_sd_rate_c;
  int<lower=-2, upper=2> theta_pow_c;

  // same things but for Theta_r
  int<lower=0> len_w7_c;
  array[Ng] int<lower=0> wg7_c;
  array[Ng] vector[len_w7_c] w7_c;
  array[Ng, len_w7_c] int<lower=1> v7_c;
  array[Ng, p_c + 1] int<lower=1> u7_c;
  array[sum(wg7_c), 3] int<lower=0> w7skel_c;
  int<lower=0> len_thet_r_c;
  array[len_thet_r_c] real<lower=0> theta_r_alpha_c;
  array[len_thet_r_c] real<lower=0> theta_r_beta_c;

  // for blocks within Theta_r that receive lkj
  array[5] int<lower=0> thetanblk_c;
  array[5] int<lower=2> thetadims_c;
  array[sum(thetanblk_c), 7] int<lower=0> thetablkse_c;
  int<lower=0> len_w8_c;
  array[Ng] int<lower=0> wg8_c;
  array[Ng] vector[len_w8_c] w8_c;
  array[Ng, len_w8_c] int<lower=1> v8_c;
  array[Ng, p_c + 1] int<lower=1> u8_c;
  array[sum(wg8_c), 3] int<lower=0> w8skel_c;
  array[Ng, p_c] int<lower=1> thetaorder_c;
  array[Ng, p_c] int<lower=1> thetarevord_c;
  
  // same things but for Psi
  int<lower=0> len_w9_c;
  array[Ng] int<lower=0> wg9_c;
  array[Ng] vector[len_w9_c] w9_c;
  array[Ng, len_w9_c] int<lower=1> v9_c;
  array[Ng, m_c + 1] int<lower=1> u9_c;
  array[sum(wg9_c), 3] int<lower=0> w9skel_c;
  int<lower=0> len_psi_sd_c;
  array[len_psi_sd_c] real<lower=0> psi_sd_shape_c;
  array[len_psi_sd_c] real<lower=0> psi_sd_rate_c;
  int<lower=-2,upper=2> psi_pow_c;
  
  // same things but for Psi_r
  int<lower=0> len_w10_c;
  array[Ng] int<lower=0> wg10_c;
  array[Ng] vector[len_w10_c] w10_c;
  array[Ng, len_w10_c] int<lower=1> v10_c;
  array[Ng, m_c + 1] int<lower=1> u10_c;
  array[sum(wg10_c), 3] int<lower=0> w10skel_c;
  array[sum(wg10_c), 3] int<lower=0> psi_r_sign_c;
  int<lower=0> len_psi_r_c;
  array[len_psi_r_c] real<lower=0> psi_r_alpha_c;
  array[len_psi_r_c] real<lower=0> psi_r_beta_c;
  int<lower=0,upper=1> fullpsi_c;

  // for blocks within Psi_r
  array[5] int<lower=0> psinblk_c;
  array[5] int<lower=3> psidims_c;
  array[sum(psinblk_c), 7] int<lower=0> psiblkse_c;
  int<lower=0> len_w11_c;
  array[Ng] int<lower=0> wg11_c;
  array[Ng] vector[len_w11_c] w11_c;
  array[Ng, len_w11_c] int<lower=1> v11_c;
  array[Ng, m_c + 1] int<lower=1> u11_c;
  array[sum(wg11_c), 3] int<lower=0> w11skel_c;
  array[Ng, m_c] int<lower=1> psiorder_c;
  array[Ng, m_c] int<lower=1> psirevord_c;
  
  // same things but for Nu
  int<lower=0> len_w13_c;
  array[Ng] int<lower=0> wg13_c;
  array[Ng] vector[len_w13_c] w13_c;
  array[Ng, len_w13_c] int<lower=1> v13_c;
  array[Ng, p_c + 1] int<lower=1> u13_c;
  array[sum(wg13_c), 3] int<lower=0> w13skel_c;
  int<lower=0> len_nu_c;
  array[len_nu_c] real nu_mn_c;
  array[len_nu_c] real<lower=0> nu_sd_c;
    
  // same things but for Alpha
  int<lower=0> len_w14_c;
  array[Ng] int<lower=0> wg14_c;
  array[Ng] vector[len_w14_c] w14_c;
  array[Ng, len_w14_c] int<lower=0> v14_c;
  array[Ng, m_c + 1] int<lower=1> u14_c;
  array[sum(wg14_c), 3] int<lower=0> w14skel_c;
  array[sum(wg14_c), 3] int<lower=0> alph_sign_c;
  int<lower=0> len_alph_c;
  array[len_alph_c] real alpha_mn_c;
  array[len_alph_c] real<lower=0> alpha_sd_c;
}
transformed data { // (re)construct skeleton matrices in Stan (not that interesting)
  array[Ng] matrix[p, m] Lambda_y_skeleton;
  array[Ng] matrix[m, m] B_skeleton;
  array[Ng] matrix[p, p] Theta_skeleton;
  array[Ng] matrix[p, p] Theta_r_skeleton;
  array[Ng] matrix[p, p] Theta_r_skeleton_f;
  array[Ng] matrix[m, m] Psi_skeleton;
  array[Ng] matrix[m, m] Psi_r_skeleton;
  array[Ng] matrix[m, m] Psi_r_skeleton_f;
  array[Ng] matrix[p, 1] Nu_skeleton;
  array[Ng] matrix[m, 1] Alpha_skeleton;
  array[Ng] matrix[sum(nlevs) - Nord, 1] Tau_skeleton;
  array[Np] vector[ord ? 0 : (p + q)] YXbarstar;
  array[Np] matrix[ord ? 0 : (p + q), ord ? 0 : (p + q)] Sstar;
  
  matrix[m, m] I = diag_matrix(rep_vector(1, m));
  
  int Ncont = p + q - Nord;
  array[max(nclus[,2]) > 1 ? max(nclus[,2]) : 0] int<lower = 0> intone;
  array[len_alph] int paidx;
  array[len_b, 2] int pbidx;
  int pridx = 1;
  int f1idx = 1;
  int f2idx = 1;
  
  array[Ng,2] int g_start1;
  array[Ng,2] int g_start4;
  array[Ng,2] int g_start5;
  array[Ng,2] int g_start7;
  array[Ng,2] int g_start9;
  array[Ng,2] int g_start10;
  array[Ng,2] int g_start13;
  array[Ng,2] int g_start14;
  array[Ng,2] int g_start15;

  array[15] int len_free;
  array[15] int pos;

  array[Ng] int matdim;
  int maxdim;
  
  for (i in 1:15) {
    len_free[i] = 0;
    pos[i] = 1;
  }

  for (g in 1:Ng) {
    Lambda_y_skeleton[g] = to_dense_matrix(p, m, w1[g], v1[g,], u1[g,]);
    B_skeleton[g] = to_dense_matrix(m, m, w4[g], v4[g,], u4[g,]);
    Theta_skeleton[g] = to_dense_matrix(p, p, w5[g], v5[g,], u5[g,]);
    Theta_r_skeleton[g] = to_dense_matrix(p, p, w7[g], v7[g,], u7[g,]);
    Theta_r_skeleton_f[g] = to_dense_matrix(p, p, w8[g], v8[g,], u8[g,]);
    Psi_skeleton[g] = to_dense_matrix(m, m, w9[g], v9[g,], u9[g,]);
    Psi_r_skeleton[g] = to_dense_matrix(m, m, w10[g], v10[g,], u10[g,]);
    Psi_r_skeleton_f[g] = to_dense_matrix(m, m, w11[g], v11[g,], u11[g,]);
    if (!use_cov) {
      Nu_skeleton[g] = to_dense_matrix((p + q), 1, w13[g], v13[g,], u13[g,]);
      Alpha_skeleton[g] = to_dense_matrix((m + n), 1, w14[g], v14[g,], u14[g,]);
    }
    Tau_skeleton[g] = to_dense_matrix(sum(nlevs) - Nord, 1, w15[g], v15[g,], u15[g,]);

    
    // count free elements in Lambda_y_skeleton
    g_start1[g,1] = len_free[1] + 1;
    g_start1[g,2] = pos[1];
    for (i in 1:p) {
      for (j in 1:m) {
        if (is_inf(Lambda_y_skeleton[g,i,j])) {
	  if (w1skel[pos[1],2] == 0 || w1skel[pos[1],3] == 1) len_free[1] += 1;
	  pos[1] += 1;
        }
      }
    }

    // same thing but for B_skeleton
    g_start4[g,1] = len_free[4] + 1;
    g_start4[g,2] = pos[4];
    for (i in 1:m) {
      for (j in 1:m) {
	if (is_inf(B_skeleton[g,i,j])) {
	  if (w4skel[pos[4],2] == 0 || w4skel[pos[4],3] == 1) len_free[4] += 1;
	  pos[4] += 1;
	}
      }
    }
    
    // same thing but for Theta_skeleton
    g_start5[g,1] = len_free[5] + 1;
    g_start5[g,2] = pos[5];
    for (i in 1:p) {
      if (is_inf(Theta_skeleton[g,i,i])) {
	if (w5skel[pos[5],2] == 0 || w5skel[pos[5],3] == 1) len_free[5] += 1;
	pos[5] += 1;
      }
    }

    // same thing but for Theta_r_skeleton
    g_start7[g,1] = len_free[7] + 1;
    g_start7[g,2] = pos[7];
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
	if (is_inf(Theta_r_skeleton[g,j,i])) {
	  if (w7skel[pos[7],2] == 0 || w7skel[pos[7],3] == 1) len_free[7] += 1;
	  pos[7] += 1;
	}
	if (is_inf(Theta_r_skeleton_f[g,j,i])) {
	  if (w8skel[pos[8],2] == 0 || w8skel[pos[8],3] == 1) len_free[8] += 1;
	  pos[8] += 1;
	}	
      }
    }

    // same thing but for Psi_skeleton
    g_start9[g,1] = len_free[9] + 1;
    g_start9[g,2] = pos[9];
    for (i in 1:m) {
      if (is_inf(Psi_skeleton[g,i,i])) {
	if (w9skel[pos[9],2] == 0 || w9skel[pos[9],3] == 1) len_free[9] += 1;
	pos[9] += 1;
      }
    }

    // same thing but for Psi_r_skeleton
    g_start10[g,1] = len_free[10] + 1;
    g_start10[g,2] = pos[10];
    for (i in 1:(m-1)) {
      for (j in (i+1):m) {
	if (is_inf(Psi_r_skeleton[g,j,i])) {
	  if (w10skel[pos[10],2] == 0 || w10skel[pos[10],3] == 1) len_free[10] += 1;
	  pos[10] += 1;
	}
	if (is_inf(Psi_r_skeleton_f[g,j,i])) {
	  if (w11skel[pos[11],2] == 0 || w11skel[pos[11],3] == 1) len_free[11] += 1;
	  pos[11] += 1;
	}
      }
    }

    if (!use_cov) {
      // same thing but for Nu_skeleton
      // pos = len_free13 + 1;
      g_start13[g,1] = len_free[13] + 1;
      g_start13[g,2] = pos[13];
      for (i in 1:(p+q)) {
	if (is_inf(Nu_skeleton[g,i,1])) {
	  if (w13skel[pos[13],2] == 0 || w13skel[pos[13],3] == 1) len_free[13] += 1;
	  pos[13] += 1;
	}
      }

      // same thing but for Alpha_skeleton
      g_start14[g,1] = len_free[14] + 1;
      g_start14[g,2] = pos[14];
      for (i in 1:(m+n)) {
	if (is_inf(Alpha_skeleton[g,i,1])) {
	  if (w14skel[pos[14],2] == 0 || w14skel[pos[14],3] == 1) len_free[14] += 1;
	  pos[14] += 1;
	}
      }
    }

    // same thing but for Tau_skeleton
    g_start15[g,1] = len_free[15] + 1;
    g_start15[g,2] = pos[15];
    for (i in 1:(sum(nlevs) - Nord)) {
      if (is_inf(Tau_skeleton[g,i,1])) {
	if (w15skel[pos[15],2] == 0 || w15skel[pos[15],3] == 1) len_free[15] += 1;
	pos[15] += 1;
      }
    }
  }

  for (g in 1:Ng) {
    if (g == Ng) {
      matdim[g] = (len_free[4] - g_start4[g, 1] + 1) + (len_free[14] - g_start14[g, 1] + 1);
    } else {
      matdim[g] = (g_start4[(g + 1), 1] - g_start4[g, 1]) + (g_start14[(g + 1), 1] - g_start14[g, 1]);
    }

    // indexing of free params across rows of Alpha combined with B
    for (r in 1:m) {
      real askel = Alpha_skeleton[g, r, 1];
      if (is_inf(askel)) {
	paidx[f1idx] = pridx;
	f1idx += 1;
	pridx += 1;
      }
      for (c in 1:m) {
	real bskel = B_skeleton[g, r, c];
	if (is_inf(bskel)) {
	  // find columnwise "free" index
	  // this could be sent in as data to improve efficiency
	  int f3idx = 1;
	  for (cc in 1:c) {
	    for (rr in 1:m) {
	      if (is_inf(B_skeleton[g, rr, cc])) {
		if (cc < c || (cc == c && rr < r)) {
		  f3idx += 1;
		}
	      }
	    }
	  }
	  pbidx[f2idx, 1] = f3idx;
	  pbidx[f2idx, 2] = pridx;
	  f2idx += 1;
	  pridx += 1;
	}
      }
    }
  }
  maxdim = max(matdim);

  if (!ord && (use_suff || use_cov)) {
    // sufficient stat matrices by pattern, moved to left for missing
    for (patt in 1:Np) {
      Sstar[patt] = rep_matrix(0, p + q, p + q);
      Sstar[patt, 1:Nobs[patt], 1:Nobs[patt]] = S[patt, Obsvar[patt, 1:Nobs[patt]], Obsvar[patt, 1:Nobs[patt]]];

      for (j in 1:Nobs[patt]) {
	YXbarstar[patt,j] = YXbar[patt, Obsvar[patt,j]];
      }
    }
  }
}
parameters {
  // free elements (possibly with inequality constraints) for coefficient matrices
  vector[len_free[1]] Lambda_y_free;
  vector<lower=0>[len_free[5]] Theta_sd_free;
  array[thetanblk[1]] corr_matrix[thetadims[1]] Theta_r_mat_1;
  array[thetanblk[2]] corr_matrix[thetadims[2]] Theta_r_mat_2;
  array[thetanblk[3]] corr_matrix[thetadims[3]] Theta_r_mat_3;
  array[thetanblk[4]] corr_matrix[thetadims[4]] Theta_r_mat_4;
  array[thetanblk[5]] corr_matrix[thetadims[5]] Theta_r_mat_5;
  vector<lower=-1,upper=1>[len_free[7]] Theta_r_free; // to use beta prior
  vector[len_free[13]] Nu_free;
  vector<lower=0>[len_free[9]] Psi_sd_tmp;
  vector[len_free[15]] Tau_ufree;
}
transformed parameters {
  array[Ng] matrix[p, m] Lambda_y;
  array[Ng] matrix[p, p] Theta_sd;
  array[Ng] matrix[p, p] T_r_lower;
  array[Ng] matrix[p, p] Theta_r;
  array[Ng] matrix[p + q, 1] Nu;
  array[Ng] matrix[m, m] Psi_tmp;

  array[Ng] matrix[sum(nlevs) - Nord, 1] Tau_un;
  array[Ng] matrix[sum(nlevs) - Nord, 1] Tau;
  vector[len_free[15]] Tau_free;
  //real tau_jacobian;
  array[Ntot] matrix[Nord, 2] Tau_by_obs;
  
  vector[len_free[1]] lambda_y_primn;
  vector[len_free[13]] nu_primn;
  vector[len_free[15]] tau_primn;
  
  array[Ng] vector[p + q] Mu;
  array[Ng] matrix[p + q, p + q] Sigma;  // model covariance matrix
  array[Ng] matrix[p + q, p + q] Sigmainv_grp;  // model covariance matrix
  array[Ng] real logdetSigma_grp;
  array[Np] matrix[p + q + 1, p + q + 1] Sigmainv;  // for updating S^-1 by missing data pattern

  array[Ntot] vector[p + q] YXstar;

  for (g in 1:Ng) {
    // model matrices
    Lambda_y[g] = fill_matrix(Lambda_y_free, Lambda_y_skeleton[g], w1skel, g_start1[g,1], g_start1[g,2]);
    Theta_sd[g] = fill_matrix(Theta_sd_free, Theta_skeleton[g], w5skel, g_start5[g,1], g_start5[g,2]);
    T_r_lower[g] = fill_matrix(Theta_r_free, Theta_r_skeleton[g], w7skel, g_start7[g,1], g_start7[g,2]);
    Theta_r[g] = T_r_lower[g] + transpose(T_r_lower[g]) - diag_matrix(rep_vector(1, p));

    // NB this handles each lv separately, as opposed to in blocks
    Psi_tmp[g] = rep_matrix(0, m, m);
    if (m > 0) {
      Psi_tmp[g] = fill_matrix(pow(Psi_sd_tmp, 2), Psi_skeleton[g], w9skel, g_start9[g,1], g_start9[g,2]);
    }
    
    if (!use_cov) {
      Nu[g] = fill_matrix(Nu_free, Nu_skeleton[g], w13skel, g_start13[g,1], g_start13[g,2]);
    }
  }

  if (sum(thetanblk) > 0) {
    // we need to define a separate parameter for each dimension of correlation matrix,
    // so we need all these Theta_r_mats
    Theta_r = fill_cov(Theta_r, thetablkse, thetanblk, Theta_r_mat_1, Theta_r_mat_2, Theta_r_mat_3, Theta_r_mat_4, Theta_r_mat_5, thetaorder, thetarevord);
  }
  
  
  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA3#v=onepage&q=LISREL&f=false
  for (g in 1:Ng) {
    if (!use_cov) {
      Mu[g] = to_vector(Nu[g]);
    } else if(has_data) {
      Mu[g] = YXbar[g]; // doesn't enter in likelihood, just for lppd + loo
    }
      
    if (p > 0) {
      Sigma[g, 1:p, 1:p] = quad_form_sym(Theta_r[g], Theta_sd[g]);
      if (m > 0) {
	Sigma[g, 1:p, 1:p] += quad_form_sym(Psi_tmp[g], transpose(Lambda_y[g]));
      }
    }
  }
  
  // obtain ordered thresholds; NB untouched for two-level models
  if (ord) {
    int opos = 1;
    int ofreepos = 1;
    //tau_jacobian = 0;
    for (g in 1:Ng) {
      int vecpos = 1;
      Tau_un[g] = fill_matrix(Tau_ufree, Tau_skeleton[g], w15skel, g_start15[g,1], g_start15[g,2]);
      for (i in 1:Nord) {
	for (j in 1:(nlevs[i] - 1)) {
	  real rc = Tau_skeleton[g, vecpos, 1];
	  int eq = w15skel[opos, 1];
	  int wig = w15skel[opos, 3];

	  if (is_inf(rc)) {
	    if (eq == 0 || wig == 1) {
	      if (j == 1) {
		Tau[g, vecpos, 1] = Tau_un[g, vecpos, 1];
	      } else {
		Tau[g, vecpos, 1] = Tau[g, (vecpos - 1), 1] + exp(Tau_un[g, vecpos, 1]);
	      }

	      Tau_free[ofreepos] = Tau[g, vecpos, 1];	      
	      // this is used if a prior goes on Tau_free, instead of Tau_ufree:
	      //if (j > 1) {
	      //  tau_jacobian += Tau_un[g, vecpos, 1]; // see https://mc-stan.org/docs/2_24/reference-manual/ordered-vector.html
	      // }
	      ofreepos += 1;
	    } else if (eq == 1) {
	      int eqent = w15skel[opos, 2];
	      Tau[g, vecpos, 1] = Tau_free[eqent];
	    }	    
	    opos += 1;
	  } else {
	    // fixed value
	    Tau[g, vecpos, 1] = Tau_un[g, vecpos, 1];
	  }	  
	  vecpos += 1;
	}
      }
    }
  }

  if (Nord > 0) {
    for (mm in 1:Np) {
      int r1 = startrow[mm];
      int r2 = endrow[mm];
      for (i in r1:r2) {
	for (j in 1:Nord) {
	  int tmpobs = YXo[i,j];
	  int vecpos = tmpobs - 1;
	  if (j > 1) {
	    vecpos += sum(nlevs[1:(j - 1)]) - j + 1;
	  }

	  if (tmpobs == 1) {
	    Tau_by_obs[i,j,1] = -30;
	    Tau_by_obs[i,j,2] = Tau[grpnum[mm], vecpos + 1, 1];
	  } else if (tmpobs == nlevs[j]) {
	    Tau_by_obs[i,j,1] = Tau[grpnum[mm], vecpos, 1];
	    Tau_by_obs[i,j,2] = 30;
	  } else {
	    Tau_by_obs[i,j,1] = Tau[grpnum[mm], vecpos, 1];
	    Tau_by_obs[i,j,2] = Tau[grpnum[mm], vecpos + 1, 1];
	  }
	}
      }
    }
  }
  
  // prior vectors
  if (wigind) {
    lambda_y_primn = fill_prior(Lambda_y_free, lambda_y_mn, w1skel);
    nu_primn = fill_prior(Nu_free, nu_mn, w13skel);
    tau_primn = fill_prior(Tau_ufree, tau_mn, w15skel);
  } else {
    lambda_y_primn = to_vector(lambda_y_mn);
    nu_primn = to_vector(nu_mn);
    tau_primn = to_vector(tau_mn);
  }

  if (Ncont > 0) {
    for (patt in 1:Np) {
      for (i in startrow[patt]:endrow[patt]) {
	for (j in 1:Ncont) {
	  YXstar[i, contidx[j]] = YX[i,j];
	}
      }
    }
  }

  // move observations to the left
  if (missing) {
    for (patt in 1:Np) {
      for (i in startrow[patt]:endrow[patt]) {
	for (j in 1:Nobs[patt]) {
	  YXstar[i,j] = YXstar[i, Obsvar[patt,j]];
	}
      }
    }
  }

  // for computing mvn with sufficient stats
  if (use_suff) {
    int blkcounter = 0;
    for (g in 1:Ng) {
      Sigmainv_grp[g] = rep_matrix(0, p + q, p + q);
      logdetSigma_grp[g] = 0;
      for (bb in 1:measnblk[g]) {
	int blkstart = measblkse[blkcounter + bb, 1];
	int blkend = measblkse[blkcounter + bb, 2];
      
	Sigmainv_grp[g, measorder[g, blkstart:blkend], measorder[g, blkstart:blkend]] = inverse_spd(Sigma[g, measorder[g, blkstart:blkend], measorder[g, blkstart:blkend]]);
	logdetSigma_grp[g] += log_determinant(Sigma[g, measorder[g, blkstart:blkend], measorder[g, blkstart:blkend]]);
      }
      blkcounter += measnblk[g];
    }
    for (patt in 1:Np) {    
      Sigmainv[patt, 1:(Nobs[patt] + 1), 1:(Nobs[patt] + 1)] = sig_inv_update(Sigmainv_grp[grpnum[patt]], Obsvar[patt,], Nobs[patt], p + q, logdetSigma_grp[grpnum[patt]]);
    }
  }
}
model { // N.B.: things declared in the model block do not get saved in the output, which is okay here

  /* transformed sd parameters for priors */
  vector[len_free[5]] Theta_pri;
  vector[len_free[9]] Psi_pri;

  /* log-likelihood */
  if (use_cov && !pri_only) {
    int blkcounter = 0;
    for (g in 1:Ng) {
      for (bb in 1:measnblk[g]) {
	int blkstart = measblkse[blkcounter + bb, 1];
	int blkend = measblkse[blkcounter + bb, 2];
	target += wishart_lpdf((N[g] - 1) * Sstar[g, measorder[g, blkstart:blkend], measorder[g, blkstart:blkend]] | N[g] - 1, Sigma[g, measorder[g, blkstart:blkend], measorder[g, blkstart:blkend]]);
      }
      blkcounter += measnblk[g];
      if (Nx[g] > 0) {
	array[Nx[g]] int xvars = Xdatvar[g, 1:Nx[g]];
	target += -wishart_lpdf((N[g] - 1) * Sstar[g, xvars, xvars] | N[g] - 1, Sigma[g, xvars, xvars]);
      }
    }
  } else if (has_data && !pri_only) {
    array[p + q] int obsidx;
    array[p + q] int xidx;
    array[p + q] int xdatidx;
    int grpidx;
    int r1;
    int r2;
    int blkcounter = 0;
        
    for (mm in 1:Np) {
      obsidx = Obsvar[mm,];
      xidx = Xvar[mm,];
      xdatidx = Xdatvar[mm,];
      grpidx = grpnum[mm];
      r1 = startrow[mm];
      r2 = endrow[mm];

      if (ord) {
	for (i in r1:r2) {
	  for (qq in 1:m) {
	    vector[ngh] marglik = log(ghwt);
	    for (gg in 1:ngh) {
	      marglik[gg] += cond_density(ghnode[gg], // * sqrt(Psi_tmp[grpidx, qq, qq]), // + Alpha[grpidx, m, 1],
					  Mu[grpidx, obsidx[1:Nobs[mm]]],
					  diagonal(Theta_sd[grpidx, obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]]),
					  Lambda_y[grpidx, obsidx[1:Nobs[mm]], qq],
					  Tau_by_obs[i], Nordobs[mm]);
	    }
	    target += log_sum_exp(marglik);
	  }
	}
      } else if (!missing) {
	for (bb in 1:measnblk[mm]) {
	  int blkstart = measblkse[blkcounter + bb, 1];
	  int blkend = measblkse[blkcounter + bb, 2];
	  if (!use_suff) {
	    target += multi_normal_lpdf(YXstar[r1:r2, measorder[mm, blkstart:blkend]] | Mu[grpidx, measorder[mm, blkstart:blkend]], Sigma[mm, measorder[mm, blkstart:blkend], measorder[mm, blkstart:blkend]]);
	  } else {
	    matrix[blkend - blkstart + 2, blkend - blkstart + 2] siginvblk = rep_matrix(0, blkend - blkstart + 2, blkend - blkstart + 2);
	    siginvblk[1:(blkend - blkstart + 1), 1:(blkend - blkstart + 1)] = Sigmainv[mm, measorder[mm, blkstart:blkend], measorder[mm, blkstart:blkend]];
	    target += multi_normal_suff(YXbarstar[mm, measorder[mm, blkstart:blkend]], Sstar[mm, measorder[mm, blkstart:blkend], measorder[mm, blkstart:blkend]], Mu[grpidx, measorder[mm, blkstart:blkend]], siginvblk, r2 - r1 + 1);
	  }
	}
	blkcounter += measnblk[mm];
	if (use_suff) {
	  // needed because we ignored the log-determinant above
	  target += -.5 * (r2 - r1 + 1) * Sigmainv[mm, Nobs[mm] + 1, Nobs[mm] + 1];
	  if (Nx[mm] > 0) {
	    target += -multi_normal_suff(YXbarstar[mm, xdatidx[1:Nx[mm]]], Sstar[mm, xdatidx[1:Nx[mm]], xdatidx[1:Nx[mm]]], Mu[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigmainv[grpidx], xidx, Nx[mm], p + q, logdetSigma_grp[grpidx]), r2 - r1 + 1);
	  }
	} else {
	  if (Nx[mm] > 0) {
	    target += -multi_normal_lpdf(YXstar[r1:r2,xdatidx[1:Nx[mm]]] | Mu[grpidx, xidx[1:Nx[mm]]], Sigma[grpidx, xidx[1:Nx[mm]], xidx[1:Nx[mm]]]);
	  }
	}
      } else {
	if (!use_suff) {
	  target += multi_normal_lpdf(YXstar[r1:r2,1:Nobs[mm]] | Mu[grpidx, obsidx[1:Nobs[mm]]], Sigma[grpidx, obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]]);

	  if (Nx[mm] > 0) {
	    target += -multi_normal_lpdf(YXstar[r1:r2,xdatidx[1:Nx[mm]]] | Mu[grpidx, xidx[1:Nx[mm]]], Sigma[grpidx, xidx[1:Nx[mm]], xidx[1:Nx[mm]]]);
	  }
	} else {
	  // sufficient stats
	  target += multi_normal_suff(YXbarstar[mm, 1:Nobs[mm]], Sstar[mm, 1:Nobs[mm], 1:Nobs[mm]], Mu[grpidx, obsidx[1:Nobs[mm]]], Sigmainv[mm, 1:(Nobs[mm] + 1), 1:(Nobs[mm] + 1)], r2 - r1 + 1);
      
	  if (Nx[mm] > 0) {
	    target += -multi_normal_suff(YXbarstar[mm, xdatidx[1:Nx[mm]]], Sstar[mm, xdatidx[1:Nx[mm]], xdatidx[1:Nx[mm]]], Mu[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigmainv[grpidx], xidx, Nx[mm], p + q, logdetSigma_grp[grpidx]), r2 - r1 + 1);
	  }
	}
      }
    }
  }
  
  /* prior densities in log-units */
  target += normal_lpdf(Lambda_y_free | lambda_y_primn, lambda_y_sd);
  target += normal_lpdf(Nu_free       | nu_primn, nu_sd);
  target += normal_lpdf(Tau_ufree      | tau_primn, tau_sd);

  /* transform sd parameters to var or prec, depending on
     what the user wants. */
  Theta_pri = Theta_sd_free;
  if (len_free[5] > 0 && theta_pow != 1) {
    for (i in 1:len_free[5]) {
      Theta_pri[i] = Theta_sd_free[i]^(theta_pow);
      target += log(abs(theta_pow)) + (theta_pow - 1)*log(Theta_sd_free[i]);
    }
  }
  Psi_pri = Psi_sd_tmp;
  if (len_free[9] > 0 && psi_pow != 1) {
    for (i in 1:len_free[9]) {
      Psi_pri[i] = Psi_sd_tmp[i]^(psi_pow);
      target += log(abs(psi_pow)) + (psi_pow - 1)*log(Psi_sd_tmp[i]);
    }
  }
  
  target += gamma_lpdf(Theta_pri | theta_sd_shape, theta_sd_rate);
  target += gamma_lpdf(Psi_pri | psi_sd_shape, psi_sd_rate);

  target += beta_lpdf(.5 * (1 + Theta_r_free) | theta_r_alpha, theta_r_beta) + log(.5) * len_free[7]; // the latter term is the jacobian moving from (-1,1) to (0,1), because beta_lpdf is defined on (0,1)
  if (sum(thetanblk) > 0) {
    for (k in 1:sum(thetanblk)) {
      int blkidx = thetablkse[k, 6];
      int arrayidx = thetablkse[k, 5];

      if (arrayidx == 1) {
	target += lkj_corr_lpdf(Theta_r_mat_1[blkidx] | thetablkse[k,7]);
      } else if (arrayidx == 2) {
	target += lkj_corr_lpdf(Theta_r_mat_2[blkidx] | thetablkse[k,7]);
      } else if (arrayidx == 3) {
	target += lkj_corr_lpdf(Theta_r_mat_3[blkidx] | thetablkse[k,7]);	
      } else if (arrayidx == 4) {
	target += lkj_corr_lpdf(Theta_r_mat_4[blkidx] | thetablkse[k,7]);
      } else {
	target += lkj_corr_lpdf(Theta_r_mat_5[blkidx] | thetablkse[k,7]);
      }      
    }
  }
}
generated quantities { // these matrices are saved in the output but do not figure into the likelihood
  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA34#v=onepage&q=LISREL&f=false

  // parameters to be obtained via Gibbs steps; all moved from earlier blocks here
  vector[len_free[4]] B_free = rep_vector(1, len_free[4]);
  vector<lower=0>[len_free[9]] Psi_sd_free = rep_vector(1, len_free[9]);
  vector[len_free[14]] Alpha_free = rep_vector(0, len_free[14]);
  array[Ng] matrix[m, m] B;
  array[Ng] matrix[m + n, 1] Alpha;
  array[Ng] matrix[m, m] Psi;
  array[Ng] matrix[m, m] PS;
  array[Ng] matrix[m, m] Psi_sd;
  array[Ng] matrix[m, m] Psi_r_lower;
  array[Ng] matrix[m, m] Psi_r;
  matrix[m, m] Psi_inv;
  array[Ng] matrix[p, p] Theta_sd_dum = Theta_sd;
  array[Ng] matrix[p, p] Thetmat;
  array[Ng] matrix[p, p] Thet;
  
  vector[len_free[9]] Psi_pri;
  vector[len_free[4]] b_primn;
  vector[len_free[14]] alpha_primn;

  array[Ntot] vector[m] eta;
  array[Ntot] vector[p + q] YXstar_gibbs;
  
  // intermediate computations for gibbs sampler
  array[Np] vector[maxdim] gamma0;
  array[Np] matrix[maxdim, maxdim] Omega_inv;
  
  // sign constraints and correlations
  vector[len_free[1]] ly_sign;
  vector[len_free[4]] bet_sign;
  vector[len_free[14]] al_sign;
  vector[len_free[8]] Theta_cov;
  vector[len_free[5]] Theta_var;
  vector[len_free[11]] Psi_cov;
  vector[len_free[9]] Psi_var;
  array[Ng] matrix[p, p] Sigma_full;
  array[Ng] vector[p] Mu_full;
  array[Ng] matrix[p + q, p + q] Sigmainv_full_grp;
  array[Ng] real logdetSigma_full_grp;
  array[Np] matrix[p + q + 1, p + q + 1] Sigmainv_full;

  // loglik + ppp
  vector[(use_cov ? Ng : Ntot)] log_lik; // for loo, etc
  vector[(use_cov ? Ng : Ntot)] log_lik_sat; // for ppp

  array[Ntot] vector[p + q] YXstar_rep; // artificial data
  vector[(use_cov ? Ng : Ntot)] log_lik_rep; // for loo, etc
  vector[(use_cov ? Ng : Ntot)] log_lik_rep_sat; // for ppp
  array[Ng] matrix[p + q, p + q + 1] satout;
  array[Ng] matrix[p + q, p + q + 1] satrep_out;
  array[Ng] vector[p + q] Mu_sat;
  array[Ng] matrix[p + q, p + q] Sigma_sat;
  array[Ng] matrix[p + q, p + q] Sigma_sat_inv_grp;
  array[Ng] real logdetS_sat_grp;
  array[Np] matrix[p + q + 1, p + q + 1] Sigma_sat_inv;
  array[Ng] vector[p + q] Mu_rep_sat;
  array[Ng] matrix[p + q, p + q] Sigma_rep_sat;
  array[Ng] matrix[p + q, p + q] Sigma_rep_sat_inv_grp;
  array[Np] matrix[p + q + 1, p + q + 1] Sigma_rep_sat_inv;
  array[Ng] real logdetS_rep_sat_grp;
  matrix[p + q, p + q] zmat;
  vector[Ng] log_lik_x_rep;
  array[Ng] matrix[p, m] Lambda;
  array[Ng] matrix[m, 1] alpha_prior;
  array[Ng] matrix[m, 1] alpha_prior_prec;
  array[Ng] matrix[m, m] b_prior;
  array[Ng] matrix[m, m] b_prior_prec;
  array[Ng] matrix[m, m] Psi_prior_shape;
  array[Ng] matrix[m, m] Psi_prior_rate;
  real<lower=0, upper=1> ppp;
    
  // Begin with Gibbs sampler of structural model
  // build prior vector/matrix for structural parameters and
  // fill structural matrices with initial values
  if (wigind) {
    b_primn = fill_prior(B_free, b_mn, w4skel);
    alpha_primn = fill_prior(Alpha_free, alpha_mn, w14skel);
  } else {
    b_primn = to_vector(b_mn);
    alpha_primn = to_vector(alpha_mn);
  }
  ly_sign = sign_constrain_load(Lambda_y_free, len_free[1], lam_y_sign);

  for (g in 1:Ng) {
    alpha_prior[g] = fill_matrix(alpha_primn, Alpha_skeleton[g], w14skel, g_start14[g,1], g_start14[g,2]);
    alpha_prior_prec[g] = fill_matrix(pow(to_vector(alpha_sd), -2), Alpha_skeleton[g], w14skel, g_start14[g,1], g_start14[g,2]);
    b_prior[g] = fill_matrix(b_primn, B_skeleton[g], w4skel, g_start4[g,1], g_start4[g,2]);
    b_prior_prec[g] = fill_matrix(pow(to_vector(b_sd), -2), B_skeleton[g], w4skel, g_start4[g,1], g_start4[g,2]);
    Psi_prior_shape[g] = fill_matrix(to_vector(psi_sd_shape), Psi_skeleton[g], w9skel, g_start9[g,1], g_start9[g,2]);
    Psi_prior_rate[g] = fill_matrix(to_vector(psi_sd_rate), Psi_skeleton[g], w9skel, g_start9[g,1], g_start9[g,2]);

    // around here, rstan line numbers are off by about 135
    Lambda[g] = fill_matrix(ly_sign, Lambda_y_skeleton[g], w1skel, g_start1[g,1], g_start1[g,2]);
    B[g] = fill_matrix(B_free, B_skeleton[g], w4skel, g_start4[g,1], g_start4[g,2]);
    Alpha[g] = fill_matrix(Alpha_free, Alpha_skeleton[g], w14skel, g_start14[g,1], g_start14[g,2]);
    Psi_r_lower[g] = fill_matrix(rep_vector(0, len_free[10]), Psi_r_skeleton[g], w10skel, g_start10[g,1], g_start10[g,2]);
    Psi_r[g] = Psi_r_lower[g] + transpose(Psi_r_lower[g]) - diag_matrix(rep_vector(1, m));
    Psi_sd[g] = fill_matrix(Psi_sd_free, Psi_skeleton[g], w9skel, g_start9[g,1], g_start9[g,2]);
    Psi[g] = quad_form_sym(Psi_r[g], Psi_sd[g]);
  }

  YXstar_gibbs = YXstar;
  for (i in 1:Ntot) {
    eta[i] = rep_vector(0, m);
  }

  // arrange prior info
  for (mm in 1:Np) {
    int pidx = 1;
    int g = grpnum[mm];

    gamma0[mm] = rep_vector(0, maxdim);
    Omega_inv[mm] = diag_matrix(gamma0[mm]);

    if (Ndum[mm] > 0) {
      for (j in 1:Ndum[mm]) {
	Theta_sd_dum[g, dum_ov_idx[mm, j], dum_ov_idx[mm, j]] = pow(.0001, .5);
      }
    }
    
    for (r in 1:m) {
      real askel = Alpha_skeleton[g, r, 1];
      if (is_inf(askel)) {
	gamma0[mm, pidx] = alpha_prior[g, r, 1];
	Omega_inv[mm, pidx, pidx] = alpha_prior_prec[g, r, 1];
	pidx += 1;
      }
      for (c in 1:m) {
	real bskel = B_skeleton[g, r, c];
	if (is_inf(bskel)) {
	  gamma0[mm, pidx] = b_prior[g, r, c];
	  Omega_inv[mm, pidx, pidx] = b_prior_prec[g, r, c];
	  pidx += 1;
	}
      }
    }
  }

  for (g in 1:Ng) {
    matrix[p, p] Thtmp = fill_matrix(Theta_r_free, Theta_r_skeleton[g], w7skel, g_start7[g,1], g_start7[g,2]);
    Thetmat[g] = Thtmp + transpose(Thtmp) - diag_matrix(rep_vector(1, p));
  }

  if (sum(thetanblk) > 0) {
    Thetmat = fill_cov(Thetmat, thetablkse, thetanblk, Theta_r_mat_1, Theta_r_mat_2, Theta_r_mat_3, Theta_r_mat_4, Theta_r_mat_5, thetaorder, thetarevord);
  }  
  
  for (i in 1:ngibbs) {
    for (mm in 1:Np) {
      array[p + q] int obsidx = Obsvar[mm,];
      array[p + q] int xidx;
      array[p + q] int xdatidx;
      matrix[m, m] IBinv;
      matrix[m, p] Lamt_Thet_inv;
      matrix[m, m] Psi0_inv;
      matrix[m, m] D;
      matrix[m, m] Dchol;
      vector[m] d;
      int r1 = startrow[mm];
      int r2 = endrow[mm];
      int g = grpnum[mm];
      vector[matdim[g]] params;
      matrix[matdim[g], matdim[g]] FVF = rep_matrix(0, matdim[g], matdim[g]);
      vector[matdim[g]] FVz = rep_vector(0, matdim[g]);
      matrix[matdim[g], matdim[g]] Dinv;
      int pidx = 1;

      IBinv = inverse(I - B[g]);
      if (Ndum_x[mm] > 0) {
	IBinv[dum_lv_x_idx[mm, 1:Ndum_x[mm]], dum_lv_x_idx[mm, 1:Ndum_x[mm]]] = rep_matrix(0, Ndum_x[mm], Ndum_x[mm]);
	for (j in 1:Ndum_x[mm]) {
	  IBinv[dum_lv_x_idx[mm, j], dum_lv_x_idx[mm, j]] = 1;
	}
      }
      
      // sample lvs
      Psi0_inv = inverse_spd( quad_form_sym(Psi[g], IBinv') );
      Lamt_Thet_inv = Lambda[g]' * inverse_spd( quad_form_sym(Thetmat[g], Theta_sd_dum[g]) ); // mxp
      
      D = inverse_spd( Lamt_Thet_inv * Lambda[g] + Psi0_inv );
      // eq (20) fsr paper
      // D = quad_form_sym(S[mm, 1:p, 1:p] - quad_form_sym(Theta_r[g], Theta_sd[g]), (Lamt_Thet_inv' * inverse_spd( Lamt_Thet_inv * Lambda[g] + Psi0_inv )));

      Dchol = cholesky_decompose(D);
      d = to_vector(Psi0_inv * IBinv * Alpha[g]);

      for (ridx in r1:r2) {
	// FIXME cannot handle combinations of ordinal and continuous
	if (ord) {
	  for (j in 1:Nord) {
	    YXstar_gibbs[ridx, obsidx[j]] = trunc_normal_rng(Nu[g, obsidx[j], 1] + Lambda[g, obsidx[j]] * eta[ridx], Theta_sd_dum[g, obsidx[j], obsidx[j]], Tau_by_obs[ridx, j, 1:2]);
	  }
	}
	eta[ridx] = multi_normal_cholesky_rng(D * (d + Lamt_Thet_inv * (YXstar_gibbs[ridx] - to_vector(Nu[g]))), Dchol);

	if (Ndum[mm] > 0) {
	  eta[ridx, dum_lv_idx[mm, 1:Ndum[mm]]] = YXstar[ridx, dum_ov_idx[mm, 1:Ndum[mm]]];
	}
      }
	
      // sample alpha, beta
      pidx = 1;
      Psi_inv = inverse_spd(Psi[g]);

      // construct F_i matrix
      // FIXME does not handle equality constraints
      for (ridx in r1:r2) {
	matrix[m, matdim[g]] etamat = rep_matrix(0, m, matdim[g]);
	vector[m] z = eta[ridx];
	pidx = 1;
	
	for (r in 1:m) {
	  real askel = Alpha_skeleton[g, r, 1];
	  if (is_inf(askel)) {
	    etamat[r, pidx] = 1;
	    pidx += 1;
	  } else if (askel != 0) {
	    z[r] += -askel;
	  }
	  for (c in 1:m) {
	    real bskel = B_skeleton[g, r, c];
	    if (is_inf(bskel)) {
	      etamat[r, pidx] = eta[ridx, c];
	      pidx += 1;
	    } else if (bskel != 0) {
	      z[r] += -bskel * eta[ridx, c];
	    }
	  }
	}
	
	FVF += etamat' * Psi_inv * etamat;
	FVz += etamat' * Psi_inv * z;
      }
      
      FVF += Omega_inv[mm, matdim[g], matdim[g]];
      FVz += Omega_inv[mm, matdim[g], matdim[g]] * gamma0[mm, matdim[g]];

      Dinv = inverse_spd(FVF);

      params = multi_normal_rng(Dinv * FVz, Dinv);

      // now put parameters in free parameter vectors
      for (j in 1:len_alph) {
	Alpha_free[j] = params[paidx[j]];
      }
      for (j in 1:len_b) {
	B_free[pbidx[j, 1]] = params[pbidx[j, 2]];
      }
    }

    // fill model matrices
    for (g in 1:Ng) {
      Alpha[g] = fill_matrix(Alpha_free, Alpha_skeleton[g], w14skel, g_start14[g,1], g_start14[g,2]);
      B[g] = fill_matrix(B_free, B_skeleton[g], w4skel, g_start4[g,1], g_start4[g,2]);
    }

    // sample Psi
    for (gg in 1:Ng) {
      int r1 = 1;
      int r2 = N[1];
      matrix[m, m] residcp = rep_matrix(0, m, m);      

      if (gg > 1) {
	r1 = sum(N[1:(gg - 1)]);
	r2 = sum(N[1:gg]);
      }

      // loop over n, get cross product of residuals (eta - (alpha + B * eta))
      for (ridx in r1:r2) {
	residcp += tcrossprod(to_matrix(eta[ridx] - (to_vector(Alpha[gg]) + B[gg] * eta[ridx])));
      }

      if (sum(psinblk) > 0) {
	matrix[m,m] Psiblk = Psi[gg, psiorder[gg], psiorder[gg]];
	matrix[m,m] residord = residcp[psiorder[gg], psiorder[gg]];
	for (k in 1:sum(psinblk)) {
	  int blkgrp = psiblkse[k, 4];	

	  if (blkgrp == gg) {
	    int srow = psiblkse[k, 1];
	    int erow = psiblkse[k, 2];
	    array[erow - srow + 1] int origrows = psirevord[gg, srow:erow];
	    
	    if (erow > srow) {
	      Psiblk[srow:erow, srow:erow] = inv_wishart_rng((r2 - r1 + 1) + Psi_prior_shape[gg, origrows[1], origrows[1]], residord[srow:erow, srow:erow] + Psi_prior_rate[gg, origrows, origrows]);
	    } else {
	      // this shouldn't happen
	      // even if Psi was originally fixed, we sample in this step
	      Psiblk[srow, srow] = inv_gamma_rng(.5 * (r2 - r1 + 1) + Psi_prior_shape[gg, origrows[1], origrows[1]], .5 * residord[srow, srow] + Psi_prior_rate[gg, origrows[1], origrows[1]]);
	    }
	  }
	}
	Psi[gg] = Psiblk[psirevord[gg], psirevord[gg]];
      }
      
      if (len_free[9] > 0) {
	for (ii in 1:m) {
	  if (is_inf(Psi_skeleton[gg, ii, ii])) {
	    Psi[gg, ii, ii] = inv_gamma_rng(.5 * (r2 - r1 + 1) + Psi_prior_shape[gg, ii, ii], .5 * residcp[ii, ii] + Psi_prior_rate[gg, ii, ii]);
	  }
	}
      }
    }
  }
  // END OF GIBBS SAMPLER
  // now deal with the rest as if it is a joint model
  
  // first deal with sign constraints:
  bet_sign = B_free; //sign_constrain_reg(B_free, len_free[4], b_sign, Lambda_y_free, Lambda_y_free);
  al_sign = Alpha_free;

  // off-diagonal covariance parameter vectors, from cor/sd matrices:
  if (p > 0  && len_free[8] > 0) {
    /* iden is created so that we can re-use cor2cov, even though
       we don't need to multiply to get covariances */
    array[Ng] matrix[p, p] iden;
    for (g in 1:Ng) {
      iden[g] = diag_matrix(rep_vector(1, p));
      Thet[g] = quad_form_sym(Thetmat[g], Theta_sd[g]);
    }
    Theta_cov = cor2cov(Thet, iden, len_free[8], Theta_r_skeleton_f, w8skel, Ng);
  }
  Theta_var = Theta_sd_free .* Theta_sd_free;

  for (g in 1:Ng) {
    PS[g] = Psi[g];
  }
  if (m > 0 && len_free[11] > 0) {
    /* iden is created so that we can re-use cor2cov, even though
       we don't need to multiply to get covariances */
    array[Ng] matrix[m, m] iden;
    for (g in 1:Ng) {
      iden[g] = diag_matrix(rep_vector(1, m));
    }
    Psi_cov = cor2cov(Psi, iden, len_free[11], Psi_r_skeleton_f, w11skel, Ng);
  }
  Psi_var = fill_vector(Psi, Psi_skeleton, w9skel, len_free[9], Ng);

  for (g in 1:Ng) {
    matrix[m, m] IBinv = inverse(diag_matrix(rep_vector(1, m)) - B[g]);
    Sigma_full[g] = quad_form_sym(Psi[g], IBinv' * Lambda[g]') + quad_form_sym(Thetmat[g], Theta_sd[g]);
    Mu_full[g] = to_vector(Nu[g] + Lambda[g] * IBinv * Alpha[g]);

    Sigmainv_full_grp[g] = inverse_spd(Sigma_full[g]);
    logdetSigma_full_grp[g] = log_determinant(Sigma_full[g]);
    for (patt in 1:Np) {    
      Sigmainv_full[patt, 1:(Nobs[patt] + 1), 1:(Nobs[patt] + 1)] = sig_inv_update(Sigmainv_full_grp[grpnum[patt]], Obsvar[patt,], Nobs[patt], p + q, logdetSigma_full_grp[grpnum[patt]]);
    }
  }

  { // log-likelihood
    array[p + q] int obsidx;
    array[p + q] int xidx;
    array[p + q] int xdatidx;
    int r1;
    int r2;
    int r3;
    int r4;
    int rr1;
    int rr2;
    int grpidx;
    int clusidx;

    if (do_test && use_cov) {
      for (g in 1:Ng) {
	Sigma_rep_sat[g] = wishart_rng(N[g] - 1, Sigma_full[g]);
      }
    } else if (do_test && has_data) {

      for (mm in 1:Np) {	
	obsidx = Obsvar[mm,];
	xidx = Xvar[mm,];
	xdatidx = Xdatvar[mm,];
	grpidx = grpnum[mm];
	r1 = startrow[mm];
	r2 = endrow[mm];

	for (jj in r1:r2) {
	  YXstar_rep[jj, 1:Nobs[mm]] = multi_normal_rng(Mu_full[grpidx, obsidx[1:Nobs[mm]]], Sigma_full[grpidx, obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]]);
	}
      }

      if (missing) {
	// start values for Mu and Sigma
	for (g in 1:Ng) {
	  Mu_sat[g] = rep_vector(0, p + q);
	  Mu_rep_sat[g] = Mu_sat[g];
	  Sigma_sat[g] = diag_matrix(rep_vector(1, p + q));
	  Sigma_rep_sat[g] = Sigma_sat[g];
	}

	for (jj in 1:emiter) {
	  satout = estep(YXstar, Mu_sat, Sigma_sat, Nobs, Obsvar, startrow, endrow, grpnum, Np, Ng);
	  satrep_out = estep(YXstar_rep, Mu_rep_sat, Sigma_rep_sat, Nobs, Obsvar, startrow, endrow, grpnum, Np, Ng);

	  // M step
	  for (g in 1:Ng) {
	    Mu_sat[g] = satout[g,,1]/N[g];
	    Sigma_sat[g] = satout[g,,2:(p + q + 1)]/N[g] - Mu_sat[g] * Mu_sat[g]';
	    Mu_rep_sat[g] = satrep_out[g,,1]/N[g];
	    Sigma_rep_sat[g] = satrep_out[g,,2:(p + q + 1)]/N[g] - Mu_rep_sat[g] * Mu_rep_sat[g]';
	  }
	}
      } else {
	// complete data; Np patterns must only correspond to groups
	for (mm in 1:Np) {
	  array[3] int arr_dims = dims(YXstar);
	  matrix[endrow[mm] - startrow[mm] + 1, arr_dims[2]] YXsmat; // crossprod needs matrix
	  matrix[endrow[mm] - startrow[mm] + 1, arr_dims[2]] YXsrepmat;
	  r1 = startrow[mm];
	  r2 = endrow[mm];
	  grpidx = grpnum[mm];
	  for (jj in 1:(p + q)) {
	    Mu_sat[grpidx,jj] = mean(YXstar[r1:r2,jj]);
	    Mu_rep_sat[grpidx,jj] = mean(YXstar_rep[r1:r2,jj]);
	  }
	  for (jj in r1:r2) {
	    YXsmat[jj - r1 + 1] = (YXstar[jj] - Mu_sat[grpidx])';
	    YXsrepmat[jj - r1 + 1] = (YXstar_rep[jj] - Mu_rep_sat[grpidx])';
	  }
	  Sigma_sat[grpidx] = crossprod(YXsmat)/N[grpidx];
	  Sigma_rep_sat[grpidx] = crossprod(YXsrepmat)/N[grpidx];
	  // FIXME? Sigma_sat[grpidx] = tcrossprod(YXsmat); does not throw an error??
	}
      }

      for (g in 1:Ng) {
	Sigma_sat_inv_grp[g] = inverse_spd(Sigma_sat[g]);
	logdetS_sat_grp[g] = log_determinant(Sigma_sat[g]);
	
	Sigma_rep_sat_inv_grp[g] = inverse_spd(Sigma_rep_sat[g]);
	logdetS_rep_sat_grp[g] = log_determinant(Sigma_rep_sat[g]);
      }

      for (mm in 1:Np) {
	Sigma_sat_inv[mm, 1:(Nobs[mm] + 1), 1:(Nobs[mm] + 1)] = sig_inv_update(Sigma_sat_inv_grp[grpnum[mm]], Obsvar[mm,], Nobs[mm], p + q, logdetS_sat_grp[grpnum[mm]]);
	Sigma_rep_sat_inv[mm, 1:(Nobs[mm] + 1), 1:(Nobs[mm] + 1)] = sig_inv_update(Sigma_rep_sat_inv_grp[grpnum[mm]], Obsvar[mm,], Nobs[mm], p + q, logdetS_rep_sat_grp[grpnum[mm]]);
      }
    }

    // compute log-likelihoods
    zmat = rep_matrix(0, p + q, p + q);

    for (mm in 1:Np) {
      obsidx = Obsvar[mm,];
      xidx = Xvar[mm, 1:(p + q)];
      xdatidx = Xdatvar[mm, 1:(p + q)];
      grpidx = grpnum[mm];
      r1 = startrow[mm];
      r2 = endrow[mm];

      if (use_cov) {
	log_lik[mm] = wishart_lpdf((N[mm] - 1) * Sstar[mm] | N[mm] - 1, Sigma_full[mm]);
	if (do_test) {
	  log_lik_sat[mm] = -log_lik[mm] + wishart_lpdf((N[mm] - 1) * Sstar[mm] | N[mm] - 1, Sstar[mm]);
	  log_lik_rep[mm] = wishart_lpdf(Sigma_rep_sat[mm] | N[mm] - 1, Sigma_full[mm]);
	  log_lik_rep_sat[mm] = wishart_lpdf(Sigma_rep_sat[mm] | N[mm] - 1, pow(N[mm] - 1, -1) * Sigma_rep_sat[mm]);
	}

	if (Nx[mm] > 0) {
	  array[Nx[mm]] int xvars = xdatidx[1:Nx[mm]];
	  log_lik[mm] += -wishart_lpdf((N[mm] - 1) * Sstar[mm, xvars, xvars] | N[mm] - 1, Sigma_full[mm, xvars, xvars]);
	  if (do_test) {
	    log_lik_sat[mm] += wishart_lpdf((N[mm] - 1) * Sstar[mm, xvars, xvars] | N[mm] - 1, Sigma_full[mm, xvars, xvars]);
	    log_lik_sat[mm] += -wishart_lpdf((N[mm] - 1) * Sstar[mm, xvars, xvars] | N[mm] - 1, Sstar[mm, xvars, xvars]);
	    log_lik_rep[mm] += -wishart_lpdf(Sigma_rep_sat[mm, xvars, xvars] | N[mm] - 1, Sigma_full[mm, xvars, xvars]);
	    log_lik_rep_sat[mm] += -wishart_lpdf(Sigma_rep_sat[mm, xvars, xvars] | N[mm] - 1, pow(N[mm] - 1, -1) * Sigma_rep_sat[mm, xvars, xvars]);
	  }
	}
      } else if (has_data) {
	for (jj in r1:r2) {
	  log_lik[jj] = multi_normal_suff(YXstar[jj, 1:Nobs[mm]], zmat[1:Nobs[mm], 1:Nobs[mm]], Mu_full[grpidx, obsidx[1:Nobs[mm]]], Sigmainv_full[mm], 1);

	  if (Nx[mm] > 0) {
	    log_lik[jj] += -multi_normal_suff(YXstar[jj, xdatidx[1:Nx[mm]]], zmat[1:Nx[mm], 1:Nx[mm]], Mu_full[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigmainv_full[grpidx], xidx, Nx[mm], p + q, logdetSigma_full_grp[grpidx]), 1);
	  }
	}
      }

      // saturated and y_rep likelihoods for ppp
      if (do_test) {
	if (!use_cov) {
	  r1 = startrow[mm];
	  r2 = endrow[mm];
	  for (jj in r1:r2) {
	    log_lik_rep[jj] = multi_normal_suff(YXstar_rep[jj, 1:Nobs[mm]], zmat[1:Nobs[mm], 1:Nobs[mm]], Mu_full[grpidx, obsidx[1:Nobs[mm]]], Sigmainv_full[mm], 1);

	    log_lik_sat[jj] = multi_normal_suff(YXstar[jj, 1:Nobs[mm]], zmat[1:Nobs[mm], 1:Nobs[mm]], Mu_sat[grpidx, obsidx[1:Nobs[mm]]], Sigma_sat_inv[mm], 1);	

	    log_lik_rep_sat[jj] = multi_normal_suff(YXstar_rep[jj, 1:Nobs[mm]], zmat[1:Nobs[mm], 1:Nobs[mm]], Mu_rep_sat[grpidx, obsidx[1:Nobs[mm]]], Sigma_rep_sat_inv[mm], 1);
	    
	    // log_lik_sat, log_lik_sat_rep
	    if (Nx[mm] > 0) {
	      log_lik_rep[jj] += -multi_normal_suff(YXstar_rep[jj, xdatidx[1:Nx[mm]]], zmat[1:Nx[mm], 1:Nx[mm]], Mu_full[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigmainv_full[grpidx], xidx, Nx[mm], p + q, logdetSigma_full_grp[grpidx]), 1);
	    
	      log_lik_sat[jj] += -multi_normal_suff(YXstar[jj, xdatidx[1:Nx[mm]]], zmat[1:Nx[mm], 1:Nx[mm]], Mu_sat[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigma_sat_inv[grpidx], xidx, Nx[mm], p + q, logdetS_sat_grp[grpidx]), 1);
	      
	      log_lik_rep_sat[jj] += -multi_normal_suff(YXstar_rep[jj, xdatidx[1:Nx[mm]]], zmat[1:Nx[mm], 1:Nx[mm]], Mu_rep_sat[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigma_rep_sat_inv[grpidx], xidx, Nx[mm], p + q, logdetS_rep_sat_grp[grpidx]), 1);
	    }
	  }
	  
	  // we subtract log_lik here so that _sat always varies and does not lead to
	  // problems with rhat and neff computations
	  log_lik_sat[r1:r2] -= log_lik[r1:r2];
	}
      }
    }
    
    if (do_test) {
      ppp = step((-sum(log_lik_rep) + sum(log_lik_rep_sat)) - (sum(log_lik_sat)));
    } else {
      ppp = 0;
    }
  }
  
} // end with a completely blank line (not even whitespace)
