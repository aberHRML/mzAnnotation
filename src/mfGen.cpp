#include <iomanip>
#include <Rcpp.h>
using namespace Rcpp;

//' Molecular Formula Generator
//' @param measured_mass accurate m/z for MF generation
//' @param max numeric \code{vector} of maximum elemental composition
//' @param min numeric \code{vector} of minimum elemental composition
//' @param tolerance mmu tolerance for MF generation
//' @param charge charge to apply to MF generation
//' @param applygr \code{boolean} denoting whether to apply the 7 golden rules
//' @return A \code{list} object containing each result as a character vector (Clean MF, MF, RDB, mass, # Isotopes, Error)
//' @export
//' @examples
//' res <- mfGen(341.10894,
//'              max=c(C = 12,iC = 0,H = 22,iH = 0,N = 0,iN = 0,O = 11,iO = 0,F = 0 ,Na = 0,
//'                    Si = 0,P = 0,S = 0,Cl = 0,iCl = 0,Br = 0,iBr = 0,K = 0,iK = 0),
//'              min=c(C = 0,iC = 0,H = 0,iH = 0,N = 0,iN = 0,O = 0,iO = 0,F = 0 ,Na = 0,
//'                    Si = 0,P = 0,S = 0,Cl = 0,iCl = 0,Br = 0,iBr = 0,K = 0,iK = 0),
//'              tolerance=0.01,charge=-1,applygr=T)
// [[Rcpp::export]]
std::vector<std::vector<std::string> > mfGen (double measured_mass, std::vector<int> max, std::vector<int> min,double tolerance, double charge,bool applygr)
{
/*
 HR2.C
 V1.02

 A program to calculate elemental compositions for a given mass.
 See the file README for details.

--------------------------------------------------------------------
 Copyright (c) 2001...2005 Joerg Hau <joerg.hau(at)dplanet.ch>.

 mail: joerg.hau@dplanet.ch
 www:  http://www.mysunrise.ch/users/joerg.hau/

 *changed version by Tobias Kind (TK), 2006 , Fiehnlab,
 *added extended valencies, added implementation of
  seven golden rules of molecular formula filtering 
 

 This program is free software; you can redistribute it and/or
 modify it under the terms of version 2 of the GNU General Public
 License as published by the Free Software Foundation. See the
 file LICENSE for details.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
--------------------------------------------------------------------

 Creation:  somewhere in 1992 by JHa.
 Revision:  2001-04-18, GPL\'d, first public release (JHa)
            2001-04-21, improved help text (JHa)
            2002-06-27, added sodium (JHa)
            2002-10-09, added 15N (JHa)
            2005-02-25, added -v option; license now GPL *v2* (JHa)
            2005-02-27, optimised code in calc loop (JHa)
            2005-02-28, verified and updated atomic masses (JHa)
            2005-06-17, added GPL text when "-h" is used (JHa)
			2006-01-01, extended version for BMC Bioinformatics publication - HR2 (TK)
			2006-03-03, added element ratio checks, extended valencies, only even electrons - HR2 (TK)
			2006-09-09,	1000x-10000x speedup hand optimized hehe. - HR2 (TK)
						-->special version for CHNSOP-F-Cl-Br-Si 
			2009-05-28, David Enot introduced the concept of \'Adducts\' (nadd) for MZedDB and corrected 
						some inaccuracies in April 2008. Manfred Beckmann has now corrected the use of 
						\'nadd\' and \'charge\' by calculating \'nadd\' from \'charge\' using abscharge = abs(charge)
						(did it manually because I couldn\'t find \'abs\') to correct \'measured_mass\' and limits, 
						but keeping nadd = 0 for rdb calculation of neutral MW (MB)
			2014-10-09, Jasen Finch (jsf9@aber.ac.uk). Integrated code for direct sourcing in R using the Rcpp package.

 This is ANSI C and should compile with any C compiler; use
 something along the lines of "gcc -Wall -O3 -o hr hr.c".
 "g++ -O2 -o myhr HR2.cpp"  on a Mac OS X G5 proc
 Optimize for speed, you may gain factor 3!
 NOW compiled under Visual C++ Express (faster than GCC) in C++ mode for boolean type.


 ---------------------------------------------------------------------
 */
            struct Functions {
            static float calc_rdb(std::vector<int> cnt_1 ,std::vector<float> val_1, int nr_el_1){
            int i;
            float sum = 2.0;
            for (i=0; i < nr_el_1; i++)
            sum += val_1[i] * cnt_1[i];
            return (sum/2.0);
            };
            static float calc_mass(std::vector<int> cnt_1 ,std::vector<double> mass_1,int charge_1, int nr_el_1){
            int i;
            double sum = 0.0;
            const double electron = 0.0005484;
            for (i=0; i < nr_el_1; i++)
            sum += mass_1[i] * cnt_1[i];
            return (sum - (charge_1 * electron));
            };
            static bool calc_element_ratios(std::vector<int> cnt_1, bool element_probability){
            bool CHNOPS_ok;  
            float HC_ratio;
            float NC_ratio;
            float OC_ratio;
            float PC_ratio;
            float SC_ratio;
            
            /* added the number of isotopes to the calculation - dle*/
            float C_count = (float)cnt_1[0]+(float)cnt_1[1]; // C_count=12C+13C
            float H_count = (float)cnt_1[2]+(float)cnt_1[3];
            float N_count = (float)cnt_1[4]+(float)cnt_1[5];
            /* modif end here */
            float O_count = (float)cnt_1[6]+(float)cnt_1[18];
            float P_count = (float)cnt_1[10];
            float S_count = (float)cnt_1[11];
            
            
            /* ELEMENT RATIOS allowed
            MIN    MAX (99.99%)
            H/C	0.07	6.00
            N/C	0.00	4.00
            O/C	0.00	3.00
            P/C	0.00	2.00
            S/C	0.00	6.00
            */	
            
            // set CHNOPS_ok = true and assume all ratios are ok
            CHNOPS_ok = true;
            /*element_probability = false;	 */
            
            
            if (C_count && H_count >0)					// C and H  must have one count anyway (remove for non-organics//
{	
            HC_ratio = H_count/C_count;
            if (element_probability)
{
            if ((HC_ratio <  0.2) || (HC_ratio >  3.0)) // this is the H/C probability check ;
            CHNOPS_ok = false;
}
            else if (HC_ratio >  6.0) // this is the normal H/C ratio check - type cast from int to float is important
            CHNOPS_ok = false;
}
            
            if (N_count >0)	// if positive number of nitrogens then thes N/C ratio else just calc normal
{
            NC_ratio = N_count/C_count;
            if (element_probability)
{
            if (NC_ratio >  2.0) // this is the N/C probability check ;
            CHNOPS_ok = false;
}
            else if (NC_ratio >  4.0)
            CHNOPS_ok = false;
}	
            
            if (O_count >0)	// if positive number of O then thes O/C ratio else just calc normal
{	
            OC_ratio = O_count/C_count;
            if (element_probability)
{
            if (OC_ratio >  1.2) // this is the O/C  probability check ;
            CHNOPS_ok = false;		
}
            else if (OC_ratio >  3.0)
            CHNOPS_ok = false;
}	
            
            if (P_count >0)	// if positive number of P then thes P/C ratio else just calc normal
{	
            PC_ratio = 	P_count/C_count;
            if (element_probability)
{
            if (PC_ratio >  0.32) // this is the P/C  probability check ;
            CHNOPS_ok = false;	
            
}
            else if (PC_ratio >  6.0)
            CHNOPS_ok = false;
}	
            
            if (S_count >0)	// if positive number of S then thes S/C ratio else just calc normal
{	
            SC_ratio = 	S_count/C_count;
            if (element_probability)
{
            if (SC_ratio >  0.65) // this is the S/C  probability check ;
            CHNOPS_ok = false;	
}
            else if (SC_ratio >  2.0)
            CHNOPS_ok = false;
}	
            
            //-----------------------------------------------------------------------------	
            
            // check for multiple element ratios together with probability check 
            //if N<10, O<20, P<4, S<3 then true
            if (element_probability && (N_count > 10) && (O_count > 20) && (P_count > 4) && (S_count > 1))
            CHNOPS_ok = false;	
            
            // NOP check for multiple element ratios together with probability check
            // NOP all > 3 and (N<11, O <22, P<6 then true)
            if (element_probability && (N_count > 3) && (O_count > 3) && (P_count > 3))
{
            if (element_probability && (N_count > 11) && (O_count > 22) && (P_count > 6))
            CHNOPS_ok = false;	
}
            
            // OPS check for multiple element ratios together with probability check
            // O<14, P<3, S<3 then true
            if (element_probability && (O_count > 14) && (P_count > 3) && (S_count > 3))
            CHNOPS_ok = false;	
            
            // PSN check for multiple element ratios together with probability check
            // P<3, S<3, N<4 then true
            if (element_probability && (P_count > 3) && (S_count > 3) && (N_count >4))
            CHNOPS_ok = false;	
            
            
            // NOS check for multiple element ratios together with probability check
            // NOS all > 6 and (N<19 O<14 S<8 then true)
            if (element_probability && (N_count >6) && (O_count >6) && (S_count >6))
{
            if (element_probability && (N_count >19) && (O_count >14) && (S_count >8))
            CHNOPS_ok = false;	
}	
            
            // function return value;
            if (CHNOPS_ok == true)
            return true;
            else 
            return false;
            }        
            };
            float mass=0;  		/* calculated mass */
            double limit_lo, limit_hi;	/* mass limits */
            float rdb, lewis, rdbori;			/* Rings & double bonds */
            long i;			
            long long hit;		/* counts the hits, with long declaration, overflow after 25h with all formulas < 2000 Da
            long = FFFFFFFFh = 4,294,967,295d*/
            int niso;
            long long counter;
            bool elementcheck;
            double error;
            int abscharge=0;
            double nadd=0;
            const int nr_el = 19; 
            
            std::string sym_1[] = {"C","iC","H","iH","N","iN","O","iO","F","Na","Si","P","S","Cl","iCl","Br","iBr","K","iK"};
            std::vector<std::string> sym(sym_1,sym_1+19);
            std::string symi_1[] = {"C","iC","H","iH","N","iN","O","iO","F","Na","Si","P","S","Cl","iCl","Br","iBr","K","iK"};
            std::vector<std::string> symi(symi_1,symi_1+19);
            double masses_1[] = {12.000000000,13.0033548378,1.0078250321,2.0141017780,14.0030740052,15.0001088984,15.9949146221,17.99916,18.99840320,22.98976967,27.9769265327,30.97376151,31.97207069,34.96885271,36.965896,78.9183376,80.916344,38.9637069,40.9618259};
            std::vector<double> masses(masses_1,masses_1+19);
            float val_1[] = {+2.0,+2.0,-1.0,-1.0,+1.0,+1.0,0.0,0.0,-1.0,-1.0,+2.0,+3.0,+4.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
            std::vector<float> val(val_1,val_1+19);
            int save_1[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            std::vector<int> save(save_1,save_1+19);
            int cnt_1[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            std::vector<int> cnt(cnt_1,cnt_1+19);
            
            std::vector<std::vector<std::string> > res;
            
            /* calculate limits */
            /* correct m/z value "measured_mass" for z>1 using "abscharge"; 
            keep charge = 0 for neutral mass searches (z=0);
            additionally, keep the idea of "adducts" for "rdb", but use charge state instead -meb */
            if (charge<0) {
            abscharge = -charge;
            nadd = abscharge;
            } else if (charge>0) {
            abscharge = charge;
            nadd = abscharge;
            } else if (charge==0) {
            abscharge = 1;			/* for limits */
            nadd = 0;				/* for rdb */
            }
            
            limit_lo = (measured_mass * abscharge) - (tolerance / 1000.0);
            limit_hi = (measured_mass * abscharge) + (tolerance / 1000.0);
            
            
            hit = 0;			/* Reset counter */
            counter = 0;
            
            /* Now lets run the big big loop ... I would like to do that
            recursively but did not yet figure out how ;-) 
            TK Adds: the loop is just fine.
            */
            
            /* now comes the "COOL trick" for calculating all formulae:
            sorting the high mass elements to the outer loops, the small weights (H)
            to the inner loops;
            
            This will reduce the computational time by factor ~10-60-1000
            OLD HR: Cangrelor at 1ppm  4465 formulas found in   5866 seconds.
            NEW HR2: Cangrelor at 1ppm 4465 formulas found in     96 seconds.
            NEW2 HR2: Cangrelor at 1ppm 4465 formulas found in     60 seconds.
            NEW3 HR2: Cangrelor at 1ppm 4465 formulas found in     59 seconds.
            HR2 Fast: Cangrelor at 1ppm 4465 formulas found in     41 seconds by evaluating 2,003,436,894 formulae.
            hr2 -c "Cangrelor" -m  774.948 -t 0.77 -C 1-64 -H 1-112 -N 0-30 -O 0-80 -P 0-12 -S 0-9 -F 0-10 -L 0-10
            
            Another additional trick is to end the 2nd.. 3rd.. 4th.. xth innermost loop
            to prevent loops which are just higher and higher in mass.
            */
            
            /*dle: process new elements */
            cnt[18] = min[18] - 1;  save[18] = cnt[18]; 
            while (cnt[18]++ < max[18]) /*"iK"*/{ 
            cnt[17] = min[17] - 1;  save[17] = cnt[17]; 
            while (cnt[17]++ < max[17]) /*"K"*/{ 
            cnt[16] = min[16] - 1;  save[16] = cnt[16]; 
            while (cnt[16]++ < max[16]) /*"iBr"*/{ 
            cnt[15] = min[15] - 1;  save[15] = cnt[15]; 
            while (cnt[15]++ < max[15]) /* "Br"*/ { 
            cnt[14] = min[14] - 1;  save[14] = cnt[14]; 
            while (cnt[14]++ < max[14]) /* "iCl"*/ { 
            cnt[13] = min[13] - 1;  save[13] = cnt[13]; 
            while (cnt[13]++ < max[13]) /*"Cl"*/ { 
            cnt[12] = min[12] - 1;  save[12] = cnt[12]; 
            while (cnt[12]++ < max[12]) /*"S"*/ { 
            cnt[11] = min[11] - 1;  save[11] = cnt[11]; 
            while (cnt[11]++ < max[11]) /*"P"*/ { 
            cnt[10] = min[10] - 1;  save[10] = cnt[10]; 
            while (cnt[10]++ < max[10]) /*"Si"*/ { 
            cnt[9] = min[9] - 1;  save[9] = cnt[9]; 
            while (cnt[9]++ < max[9]) /*"Na"*/ { 
            cnt[8] = min[8] - 1;  save[8] = cnt[8]; 
            while (cnt[8]++ < max[8]) /*"F"*/{ 
            cnt[7] = min[7] - 1;  save[7] = cnt[7]; 
            while (cnt[7]++ < max[7]) /*"iO"*/ { 
            cnt[6] = min[6] - 1;  save[6] = cnt[6]; 
            while (cnt[6]++ < max[6]) /*"O"*/ { 
            cnt[5] = min[5] - 1;  save[5] = cnt[5]; 
            while (cnt[5]++ < max[5]) /*"15N"*/{ 
            cnt[4] = min[4] - 1;  save[4] = cnt[4]; 
            while (cnt[4]++ < max[4]) /*"N"*/{ 
            cnt[1] = min[1] - 1; save[1] = cnt[1]; 
            while (cnt[1]++ < max[1]) /*"13C"*/ { 
            cnt[0] = min[0] - 1; save[0] = cnt[0]; 
            while (cnt[0]++ < max[0]) /* "C"*/ { 
            cnt[3] = min[3] - 1; 	save[3] = cnt[3]; 
            while (cnt[3]++ < max[3]) /*"D"*/{ 
            cnt[2] = min[2] - 1; save[2] = cnt[2]; 
            while (cnt[2]++ < max[2]) /*"H"*/{ 
            
            mass = Functions::calc_mass(cnt,masses,charge,nr_el);
            counter++;
            
            //just for debug purposes
            //if (mass > limit_hi)  
            //printf("mass: %f\\tC: %d  H: %d  N: %d O: %d P: %d S: %d Cl: %d Br: %d\\n",mass,el[0].cnt,el[2].cnt,el[4].cnt,el[6].cnt,el[10].cnt,el[11].cnt,el[12].cnt,el[13].cnt);
            
            /* if we exceed the upper limit, we can stop the calculation
            for this particular element (JHa 20050227). <-- comment TK that will only bust the innermost while loop, which is "H"*/
            
            // break H loop 
            if (mass > limit_hi)  break;
            //************************************************************************************************************/	
            //Calculus loop with print out
            //************************************************************************************************************/	
            
            if ((mass >= limit_lo) && (mass <= limit_hi)){ /* within limits? */	
            // element check will be performed always, if variable bool element_probability is true also probabilities will be calculated
            // not an elegant implementation, but fast.
            elementcheck = Functions::calc_element_ratios(cnt,applygr);  /* pass applygr boolean by dle */
            if (elementcheck){ 
            rdbori = Functions::calc_rdb(cnt,val,nr_el);	/* get RDB */           
            rdb = rdbori + 0.5*nadd; /* dle: if nadd addcuts */
            lewis = (float)(fmod(rdb, 1)); /*calc reminder*/
            if ((rdb >= 0) && (lewis != 0.5) && (lewis !=-0.5)){ /* less than -0.5 RDB does not make sense */													/* NO(!) CH3F10NS2 exists , RDB =  -4.0   M= 282.9547*/
            hit ++;
            std::string clean_mf;
            for (i = 0; i < nr_el; i++){			/* print composition */
            if (cnt[i] > 0) {				/* but only if useful */
            clean_mf += sym[i];
            std::string count;
            std::ostringstream convert;
            convert << cnt[i];
            count = convert.str();
            clean_mf += count;
            clean_mf += ".";
            }
            }
            /* dle: print out a more explicit molecular formula for further processing and
            variable niso counts number of isotope elements in the solution */
            niso=0;
            std::string mf;
            for (i = 0; i < nr_el; i++)
{			/* print composition */
            if (cnt[i] > 0)
{				/* but only if useful */
            mf += symi[i];
            std::string count;
            std::ostringstream convert;
            convert << cnt[i];
            count = convert.str();
            mf += count;
            if (symi[i] != sym[i])
            niso=niso+cnt[i];
}
            
}
            /* dle: end of molecular print out */
            
            error = 1000.0 * ((measured_mass * abscharge) - mass);
            /* convert results to strings for output */
            std::string out_rdb;
            std::string out_mass;
            std::string out_nadd;
            std::string out_niso;
            std::string out_error;
            std::ostringstream convert_1;
            convert_1 << rdb;
            out_rdb = convert_1.str();
            std::ostringstream convert_2;
            convert_2 << std::setprecision(9) << mass;
            out_mass = convert_2.str();
            std::ostringstream convert_3;
            convert_3 << nadd;
            out_nadd = convert_3.str();
            std::ostringstream convert_4;
            convert_4 << niso;
            out_niso = convert_4.str();
            std::ostringstream convert_5;
            convert_5 << error;
            out_error = convert_5.str();
            
            std::vector<std::string> out(7);
            out[0] = clean_mf;
            out[1] = mf;
            out[2] = out_rdb;
            out[3] = out_mass;
            out[4] = out_niso;
            out[5] = out_error;
            res.push_back(out);
            }	/* end of "rdb" loop */
            
            }	// end of elementcheck loop
            
            }	/* end of "limit" loop */
            //************************************************************************************************************/
            /*
            TK: if the current mass is larger than the limit the loop can be exited.
            Each element must point to the element which is in use and before.
            This is a static implementation which can be enhanced with a pointer chain to the lower element.
            Actually now its only allowed for CHNSOP-Fl-Cl-Br-Si !!! Brute-force <> elegance :-)
            */
            } /*"H"*/
            if ((mass >= limit_lo) && (save[2] == cnt[2]-1)) break;
            } /*"D"*/
            if ((mass >= limit_lo) && (save[3] == cnt[3]-1)) break;  /* dle addons */
            } /* "C"*/
            if ((mass >= limit_lo) && (save[0] == cnt[0]-1)) break;
            } /*"13C"*/
            if ((mass >= limit_lo) && (save[1] == cnt[1]-1)) break;  /* dle addons */
            } /*"N"*/
            if ((mass >= limit_lo) && (save[4] == cnt[4]-1)) break;
            } /*"15N"*/
            if ((mass >= limit_lo) && (save[5] == cnt[5]-1)) break;  /* dle addons */
            } /*"O"*/
            if ((mass >= limit_lo) && (save[6] == cnt[6]-1)) break;
            } /*"iO"*/
            if ((mass >= limit_lo) && (save[7] == cnt[7]-1)) break;
            } /*"F"*/
            if ((mass >= limit_lo) && (save[8] == cnt[8]-1)) break;
            }  /*"Na"*/
            if ((mass >= limit_lo) && (save[9] == cnt[9]-1)) break;
            } /*"Si"*/
            if ((mass >= limit_lo) && (save[10] == cnt[10]-1)) break;
            } /*"P"*/
            if ((mass >= limit_lo) && (save[11] == cnt[11]-1)) break; 
            } /*"S"*/
            if ((mass >= limit_lo) && (save[12] == cnt[12]-1)) break;
            } /*"Cl"*/
            if ((mass >= limit_lo) && (save[13] == cnt[13]-1)) break;  /* dle addons */
            } /*"iCL"*/
            if ((mass >= limit_lo) && (save[14] == cnt[14]-1)) break;
            } /*"Br"*/
            if ((mass >= limit_lo) && (save[15] == cnt[15]-1)) break;	 /* dle addons */
            } /*"iBr"*/
            if ((mass >= limit_lo) && (save[16] == cnt[16]-1)) break;	 /* dle addons */
            } /*"K"*/
            if ((mass >= limit_lo) && (save[17] == cnt[17]-1)) break;	 /* dle addons */
            } /*"iK"*/
            /* close that giant loop thing started above */
            return res;
}
