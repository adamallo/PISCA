	/*
	 * LOH.java
	 *
	 * Modified from HKY.java by DM
	 *
	 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
	 *
	 * This file is part of BEAST.
	 * See the NOTICE file distributed with this work for additional
	 * information regarding copyright ownership and licensing.
	 *
	 * BEAST is free software; you can redistribute it and/or modify
	 * it under the terms of the GNU Lesser General Public License as
	 * published by the Free Software Foundation; either version 2
	 * of the License, or (at your option) any later version.
	 *
	 *  BEAST is distributed in the hope that it will be useful,
	 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
	 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	 *  GNU Lesser General Public License for more details.
	 *
	 * You should have received a copy of the GNU Lesser General Public
	 * License along with BEAST; if not, write to the
	 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
	 * Boston, MA  02110-1301  USA
	 */
	
package PISCA;
	
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evolution.datatype.TwoStates;
import PISCA.LOHParser;
import dr.evomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
	
	
	/**
	 * LOH model
	 *
	 * @version $Id: LOH.java $
	 *
	 * @author Rumen Kostadinov
	 * @author Diego Mallo
	 */
public class LOH extends AbstractSubstitutionModel implements dr.util.XHTMLable {
	
	
	    private Variable<Double> alphaParameter = null;
	    private Variable<Double> betaParameter = null;
	
	    private boolean updateIntermediates = true;
	    
	    private double alpha=0;
	    private double beta=0;
	
	    /**
	     * Constructor
	     */
	    public LOH(Variable alphaParameter, Variable betaParameter,FrequencyModel freqModel) {
	
	        super(LOHParser.LOH_MODEL,TwoStates.INSTANCE,freqModel);
	        this.alphaParameter = alphaParameter;
	        addVariable(alphaParameter);
	        alphaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
	        this.betaParameter = betaParameter;
	        addVariable(betaParameter);
	        betaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
	        
	        updateIntermediates = true;
	    }
	
	    /**
	     * set alpha
	     */
	    public void setAlpha(double alpha) {
	        alphaParameter.setValue(0, alpha);
	        updateMatrix = true;
	    }
	    
	    /**
	     * set beta
	     */
	    public void setBeta(double beta) {
	        betaParameter.setValue(0, beta);
	        updateMatrix = true;
	    }
	
	    /**
	     * @return alpha
	     */
	    public final double getAlpha() {
	        return alphaParameter.getValue(0);
	    }
	    
	    /**
	     * @return beta
	     */
	    public final double getBeta() {
	        return betaParameter.getValue(0);
	    }
	
	    protected void frequenciesChanged() {
	        // frequencyModel changed
	        //updateIntermediates = true; //DM: LOH model does not have frequencies
	    }

	    protected void ratesChanged() {
	        // Nothing to precalculate
	    		
	    }
	
	    private void calculateIntermediates() {
	            
	        //final double r1 = (1 / freqR) - 1;
	        //tab1A = freqA * r1;
	
	        //tab3A = freqA / freqR;
	        //tab2A = 1 - tab3A;        // (freqR-freqA)/freqR;
	
	        //final double r2 = 1 / r1; // ((1 / freqY) - 1);
	        //tab1C = freqC * r2;
	
	        //tab3C = freqC / freqY;
	        //tab2C = 1 - tab3C;       // (freqY-freqC)/freqY; assert  tab2C + tab3C == 1.0;
	
	        //tab1G = freqG * r1;
	        //tab3G = tab2A;            // 1 - tab3A; // freqG/freqR;
	        //tab2G = tab3A;            // 1 - tab3G; // (freqR-freqG)/freqR;
	
	        //tab1T = freqT * r2;
	
	        //tab3T = tab2C;            // 1 - tab3C;  // freqT/freqY;
	        //tab2T = tab3C;            // 1 - tab3T; // (freqY-freqT)/freqY; //assert tab2T + tab3T == 1.0 ;
	    	
	        updateMatrix = true;
	        updateIntermediates = false;
	    }
	
	    /**
	     * get the complete transition probability matrix for the given distance
	     *
	     * @param distance the expected number of substitutions
	     * @param matrix   an array to store the matrix
	     */
	    public void getTransitionProbabilities(double distance, double[] matrix) {
	        synchronized (this) {
	            if (updateIntermediates) {
	                calculateIntermediates();
	            }
	
	            if (updateMatrix) {
	                setupMatrix();
	            }
	        }
	        
	        //final double xx = beta * distance;
	        //final double bbR = Math.exp(xx * A_R);
	        //final double bbY = Math.exp(xx * A_Y);
	
	       //final double aa = Math.exp(xx);
	        //final double oneminusa = 1 - aa;
	
	        //final double t1Aaa = (tab1A * aa);
	        //matrix[0] = freqA + t1Aaa + (tab2A * bbR);
	
	        //matrix[1] = freqC * oneminusa;
	        //final double t1Gaa = (tab1G * aa);
	        //matrix[2] = freqG + t1Gaa - (tab3G * bbR);
	        //matrix[3] = freqT * oneminusa;
	
	        //matrix[4] = freqA * oneminusa;
	        //final double t1Caa = (tab1C * aa);
	        //matrix[5] = freqC + t1Caa + (tab2C * bbY);
	        //matrix[6] = freqG * oneminusa;
	        //final double t1Taa = (tab1T * aa);
	        //matrix[7] = freqT + t1Taa - (tab3T * bbY);
	
	        //matrix[8] = freqA + t1Aaa - (tab3A * bbR);
	        //matrix[9] = matrix[1];
	        //matrix[10] = freqG + t1Gaa + (tab2G * bbR);
	        //matrix[11] = matrix[3];
	
	        //matrix[12] = matrix[4];
	        //matrix[13] = freqC + t1Caa - (tab3C * bbY);
	        //matrix[14] = matrix[6];
	        //matrix[15] = freqT + t1Taa + (tab2T * bbY);
	        
	        // LOH Model matrix
	        // from/to       0     1
	        //         0    1-a    a
	        //         1     b    1-b
	        alpha = getAlpha();
	        beta = getBeta();  
	        matrix[1] = alpha * (1-Math.exp(-(alpha+beta)*distance)) / (alpha+beta);
	        matrix[2] = beta  * (1-Math.exp(-(alpha+beta)*distance)) / (alpha+beta);
	        matrix[0] = 1 - matrix[1];
	        matrix[3] = 1 - matrix[2];
	
	    }
	
	    /**
	     * setup substitution matrix
	     */
	    public void setupMatrix() {
	        //final double kappa = getKappa();
	        //beta = -1.0 / (2.0 * (freqR * freqY + kappa * (freqA * freqG + freqC * freqT)));
	
	        //A_R = 1.0 + freqR * (kappa - 1);
	        //A_Y = 1.0 + freqY * (kappa - 1); 
    	
	        updateMatrix = false;
	    }
	
	    protected void setupRelativeRates() {
	    }
	
	    /**
	     * This function returns the Eigen vectors.
	     *
	     * @return the array
	     */
//	    public double[][] getEigenVectors() {
//        synchronized (this) {
//            if (updateMatrix) {
//                setupEigenSystem();
//            }
//        }
//        return Evec;
//    }
//
//    /**
//     * This function returns the inverse Eigen vectors.
//     *
//     * @return the array
//     */
//    public double[][] getInverseEigenVectors() {
//        synchronized (this) {
//            if (updateMatrix) {
//                setupEigenSystem();
//            }
//        }
//        return Ievc;
//    }
//
//    /**
//     * This function returns the Eigen values.
//     */
//    public double[] getEigenValues() {
//        synchronized (this) {
//            if (updateMatrix) {
//                setupEigenSystem();
//            }
//        }
//        return Eval;
//    }
//
//    /**
//     * setup substitution matrix
//     */
//    protected void setupEigenSystem() {
//        if (!eigenInitialised)
//            initialiseEigen();
//
//        final double kappa = getKappa();
//
//        // left eigenvector #1
//        Ievc[0][0] = freqA; // or, evec[0] = pi;
//        Ievc[0][1] = freqC;
//        Ievc[0][2] = freqG;
//        Ievc[0][3] = freqT;
//
//        // left eigenvector #2
//        Ievc[1][0] =  freqA * freqY;
//        Ievc[1][1] = -freqC * freqR;
//        Ievc[1][2] =  freqG * freqY;
//        Ievc[1][3] = -freqT * freqR;
//
//        Ievc[2][1] =  1; // left eigenvectors 3 = (0,1,0,-1); 4 = (1,0,-1,0)
//        Ievc[2][3] = -1;
//
//        Ievc[3][0] =  1;
//        Ievc[3][2] = -1;
//
//        Evec[0][0] =  1; // right eigenvector 1 = (1,1,1,1)'
//        Evec[1][0] =  1;
//        Evec[2][0] =  1;
//        Evec[3][0] =  1;
//
//        // right eigenvector #2
//        Evec[0][1] =  1.0/freqR;
//        Evec[1][1] = -1.0/freqY;
//        Evec[2][1] =  1.0/freqR;
//        Evec[3][1] = -1.0/freqY;
//
//        // right eigenvector #3
//        Evec[1][2] =  freqT / freqY;
//        Evec[3][2] = -freqC / freqY;
//
//        // right eigenvector #4
//        Evec[0][3] =  freqG / freqR;
//        Evec[2][3] = -freqA / freqR;
//
//        // eigenvectors
//        beta = -1.0 / (2.0 * (freqR * freqY + kappa * (freqA * freqG + freqC * freqT)));
//
//        A_R = 1.0 + freqR * (kappa - 1);
//        A_Y = 1.0 + freqY * (kappa - 1);
//
//        Eval[1] = beta;
//        Eval[2] = beta*A_Y;
//        Eval[3] = beta*A_R;
//
//        updateMatrix = false;
//    }
//
//    /**
//     * allocate memory for the Eigen routines
//     */
//    protected void initialiseEigen() {
//
//        Eval = new double[stateCount];
//        Evec = new double[stateCount][stateCount];
//        Ievc = new double[stateCount][stateCount];
//
//        eigenInitialised = true;
//        updateMatrix = true;
//    }
	    
	    // *****************************************************************
	    // Interface Model
	    // *****************************************************************


	    protected void storeState() {
	    } // nothing to do

	    /**
	     * Restore the additional stored state
	     */
	    protected void restoreState() {
	        updateMatrix = true;
	    }

	    protected void acceptState() {
	    } // nothing to do
	
	    // **************************************************************
	    // XHTMLable IMPLEMENTATION
	    // **************************************************************
	
	    public String toXHTML() {
	        StringBuffer buffer = new StringBuffer();
	
	        buffer.append("<em>LOH Model</em> (alpha= ");
	        buffer.append(getAlpha());
	        buffer.append(", beta = ");
	        buffer.append(getBeta());
	        buffer.append(")");
	
	        return buffer.toString();
	    };
}
