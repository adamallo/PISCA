/*
 * AbsoluteCNASubstitutionmodel.java
 *
 * By Diego Mallo
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
import PISCA.AbsoluteCNASubstitutionModelParser;
import dr.evolution.datatype.DataType;
import dr.inference.model.DuplicatedParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.AbstractModel;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.math.matrixAlgebra.RobustSingularValueDecomposition;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.Property;

/**
 * Eigen-decomposition framework copied from ComplexSubstitutionModel.java :A general irreversible class for any
 * data type; allows complex eigenstructures.
 *
 * @author Diego Mallo
 * @author Marc Suchard
 */

public class AbsoluteCNASubstitutionModel extends AbstractModel implements SubstitutionModel, dr.util.XHTMLable {
    
	public static final String MODEL = "model";

    protected DataType dataType = null;

    protected FrequencyModel freqModel;
    protected double[] relativeRates;

    public double[] getRelativeRates() {
        return relativeRates;
    }

    protected double[] storedRelativeRates;

    protected int stateCount;
    protected int rateCount;

    protected boolean eigenInitialised = false;
    protected boolean updateMatrix = true;
    protected boolean storedUpdateMatrix = true;
    protected boolean wellConditioned = true;
    private boolean storedWellConditioned = true;
    private boolean isComplex = false;
    
    private double maxConditionNumber = 10000;
    private int maxIterations = 1000;
    private boolean checkConditioning = true;
    private boolean doNormalization = false; 
	//We do not normalize and make the rates relative to the total rate when only one allele is present
    //This means that the clock rate is given in alterations per fragment per allele.
    //I do not think we can normalize with the flux of out states due to the lack of equilibrium (0 is a absorbing state).
	
	private Variable<Double> gainParameter; //Gain rate in gains per fragment per allele. Typically, this will be fixed at 1, to use a clock rate. ClockRate=gainParameter+gainParameter*lossParameter. Which will give when fixing gain at 1; gainRate=clockrate/(1+lossParameter); lossRate=clockRate-gainRate.
	private Variable<Double> lossParameter; //Relative loss rate (multiplier of the gain rate)
    
    public AbsoluteCNASubstitutionModel(DataType dataType, Variable gainParameter, Variable lossParameter, FrequencyModel freqModel) {
		super(AbsoluteCNASubstitutionModelParser.AbsoluteCNA_MODEL);
		this.dataType = dataType;

	    setStateCount(dataType.getStateCount());

	    if (freqModel != null) {
	    		if (freqModel.getDataType() != dataType) {
	                throw new IllegalArgumentException("Datatypes do not match!");
	            }

	    			this.freqModel = freqModel;
	            addModel(freqModel);
	            
	    }
		
		if (gainParameter != null) {
            addVariable(gainParameter);
            if (!(gainParameter instanceof DuplicatedParameter))
                gainParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.gainParameter = gainParameter;
		}
		if (lossParameter != null) {
            addVariable(lossParameter);
            if (!(lossParameter instanceof DuplicatedParameter))
                lossParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.lossParameter = lossParameter;
		}
	    updateMatrix = true;
    }

	protected void handleModelChangedEvent(Model model, Object object, int index) { //Executed when a registered model changes
    		if (model == freqModel) //The only model that can be registered is the frequencyModel, which is not used by this substitution model
    			return;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) { //Executed when a registered variable changes
        // Rates changed
    		updateMatrix = true;
    }
    
    private void setStateCount(int stateCount) {
        eigenInitialised = false;

        this.stateCount = stateCount;
        rateCount = (stateCount - 1) * stateCount;

        relativeRates = new double[rateCount];
        storedRelativeRates = new double[rateCount];
        for (int i = 0; i < rateCount; i++) {
            relativeRates[i] = 1.0;
        }
    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        updateMatrix = storedUpdateMatrix;
        wellConditioned = storedWellConditioned;

        // To restore all this stuff just swap the pointers...
        double[] tmp1 = storedRelativeRates;
        storedRelativeRates = relativeRates;
        relativeRates = tmp1;

        tmp1 = storedEval;
        storedEval = Eval;
        Eval = tmp1;

        double[][] tmp2 = storedIevc;
        storedIevc = Ievc;
        Ievc = tmp2;

        tmp2 = storedEvec;
        storedEvec = Evec;
        Evec = tmp2;
        
        double[] tmp3 = storedEvalImag;
        storedEvalImag = EvalImag;
        EvalImag = tmp3;

    }

    protected void storeState() {

        storedUpdateMatrix = updateMatrix;
        storedWellConditioned = wellConditioned;
        
        System.arraycopy(relativeRates, 0, storedRelativeRates, 0, rateCount);

        System.arraycopy(Eval, 0, storedEval, 0, stateCount);
        System.arraycopy(EvalImag, 0, storedEvalImag, 0, stateCount);
        
        for (int i = 0; i < stateCount; i++) {
            System.arraycopy(Ievc[i], 0, storedIevc[i], 0, stateCount);
            System.arraycopy(Evec[i], 0, storedEvec[i], 0, stateCount);
        }

    }
	
    /**
     * setup substitution matrix
     */
    public void setupMatrix() {

        if (!eigenInitialised)
        {
            initialiseEigen();
        }

        int i, j, k = 0;

        storeIntoAmat();
        
        makeValid(amat, stateCount);
        
        // compute eigenvalues and eigenvectors for reversible matrix
        
        /*elmhes(amat, ordr, stateCount); //upper hessemberg (returns ordr)
        eltran(amat, Evec, ordr, stateCount);  //Accumulates elmhes similarity transformations (returns Evec)
        hqr2(stateCount, 1, stateCount, amat, Evec, Eval, evali); //Find Evec and Evac of an upper hessemberg by QR method
        luinverse(Evec, Ievc, stateCount); // Generates Ivec, the inverse of Evec
         */
        // compute eigenvalues and eigenvectors
        RobustEigenDecomposition eigenDecomp;
        try {
            eigenDecomp = new RobustEigenDecomposition(new DenseDoubleMatrix2D(amat), maxIterations);
        } catch (ArithmeticException ae) {
            System.err.println(ae.getMessage());
            wellConditioned = false;
            System.err.println("amat = \n" + new Matrix(amat));
            return;
        }

        DoubleMatrix2D eigenV = eigenDecomp.getV();
        DoubleMatrix1D eigenVReal = eigenDecomp.getRealEigenvalues();
        DoubleMatrix1D eigenVImag = eigenDecomp.getImagEigenvalues();
        DoubleMatrix2D eigenVInv;

        // A better (?) approach to checking diagonalizability comes from:
        //
        // J. Gentle (2007) Matrix Algebra
        //
        // Diagonalizbility Theorem: A matrix A is (complex) diagonalizable iff all distinct eigenvalues \lambda_l
        // with algebraic multiplicity m_l are semi-simple, i.e.
        //
        //          rank(A - \lambda_l I) = n - m_l
        //
        // Equivalently (?), eigenV must be non-singular.
        //
        // SVD is needed to numerically approximate the rank of a matrix, so we can check Algrebra.rank()
        // or Algebra.cond() with almost equal amounts of work.  I don't know which is more reliable. -- MAS

        if (checkConditioning) {
            RobustSingularValueDecomposition svd;
            try {
                svd = new RobustSingularValueDecomposition(eigenV, maxIterations);
            } catch (ArithmeticException ae) {
                System.err.println(ae.getMessage());
                wellConditioned = false;
                return;
            }
            if (svd.cond() > maxConditionNumber) {
                wellConditioned = false;
                return;
            }
        }

        try {
            eigenVInv = alegbra.inverse(eigenV);
        } catch (IllegalArgumentException e) {
            wellConditioned = false;
            return;
        }

        Ievc = eigenVInv.toArray();
        Evec = eigenV.toArray();
        Eval = eigenVReal.toArray();
        EvalImag = eigenVImag.toArray();

        // Check for valid decomposition
        for (i = 0; i < stateCount; i++) {
            if (Double.isNaN(Eval[i]) || Double.isNaN(EvalImag[i]) ||
                    Double.isInfinite(Eval[i]) || Double.isInfinite(EvalImag[i])) {
                wellConditioned = false;
                return;
            } else if (Math.abs(Eval[i]) < 1e-10) {
                Eval[i] = 0.0;
            }
        }

        updateMatrix = false;
        wellConditioned = true;
        // compute normalization and rescale eigenvalues

        //computeStationaryDistribution();

        if (doNormalization) {
            double subst = 0.0;

            for (i = 0; i < stateCount; i++)
                subst += -amat[i][i]; //* stationaryDistribution[i];

//        normalization = subst;

            for (i = 0; i < stateCount; i++) {
                Eval[i] /= subst/stateCount;//Unweighted mean
                EvalImag[i] /= subst/stateCount;
            }
        }
    }

    //store the infinitesimal rates in the vector to a matrix called amat
    public void storeIntoAmat(){
    	

        double tg = gainParameter.getValue(0);
        double td = lossParameter.getValue(0);
        
        //We cannot normalize the matrix using the flux of states, since we do not have an equilibrium or equilibrium frequencies.
        //However, we can normalize the rates making them relative. This way we do have a general mutation rate per fragment per allele, composed by the two relative rates of gain and loss
        //This way we can still use clock models
        
        double totalRate=tg+(tg*td);
        double g = tg/totalRate;
        double d = 1-g;
        
        double [][] r = amat;
        int i, j = 0;
        
        for (i = 0; i < stateCount; i++) {
            for (j = 0; j < stateCount; j++) {
                r[i][j] = 0;
            }
        }
        
		//We skip the first row because it is all 0s
		//We are filling up the rates programatically. One to the left from the diagonal is n*d and right n*g
        //We leave the final row out of the loop, since it does not follow the same pattern
		for (i = 1; i < stateCount-1; i++) {
			r[i][i-1]=d*i;
			r[i][i+1]=g*i;
        }
		
		//Last row, only deletions allowed
		r[stateCount-1][stateCount-2]=d*(stateCount-1);
        

    }
    
    protected void acceptState() {
    } // nothing to do
/*

    protected void setupRelativeRates() {
    	
    }*/

    /**
     * @return the frequency model
     */
    public FrequencyModel getFrequencyModel() {
        return freqModel;
    }

    /**
     * @return the data type
     */
    public DataType getDataType() {
        return dataType;
    }

    /**
     * get the complete transition probability matrix for the given distance
     *
     * @param distance the expected number of substitutions
     * @param matrix   an array to store the matrix
     */
    public void getTransitionProbabilities(double distance, double[] matrix) {
        int i, j, k;
        double temp;

        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads - AJD
        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }
        
        if (!wellConditioned) {
            Arrays.fill(matrix, 0.0);
            return;
        }

        // implemented a pool of iexp matrices to support multiple threads
        // without creating a new matrix each call. - AJD
        double[][] iexp = popiexp();
        
        // Eigenvalues and eigenvectors of a real matrix A.
        //
        // If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal
        // and the eigenvector matrix V is orthogonal. I.e. A = V D V^t and V V^t equals
        // the identity matrix.
        //
        // If A is not symmetric, then the eigenvalue matrix D is block diagonal with
        // the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
        // lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda]. The columns
        // of V represent the eigenvectors in the sense that A*V = V*D. The matrix
        // V may be badly conditioned, or even singular, so the validity of the
        // equation A = V D V^{-1} depends on the conditioning of V.


        for (i = 0; i < stateCount; i++) {

                    if (EvalImag[i] == 0) {
                        // 1x1 block
                        temp = Math.exp(distance * Eval[i]);
                        for (j = 0; j < stateCount; j++) {
                            iexp[i][j] = Ievc[i][j] * temp;
                        }
                    } else {
                        // 2x2 conjugate block
                        // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                        // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                        int i2 = i + 1;
                        double b = EvalImag[i];
                        double expat = Math.exp(distance * Eval[i]);
                        double expatcosbt = expat * Math.cos(distance * b);
                        double expatsinbt = expat * Math.sin(distance * b);

                        for (j = 0; j < stateCount; j++) {
                            iexp[i][j] = expatcosbt * Ievc[i][j] + expatsinbt * Ievc[i2][j];
                            iexp[i2][j] = expatcosbt * Ievc[i2][j] - expatsinbt * Ievc[i][j];
                        }
                        i++; // processed two conjugate rows
                    }
        }

        int u = 0;
        for (i = 0; i < stateCount; i++) {
            for (j = 0; j < stateCount; j++) {
                temp = 0.0;
                for (k = 0; k < stateCount; k++) {
                    temp += Evec[i][k] * iexp[k][j];
                }
                if (temp < 0.0)
                    matrix[u] = minProb;
                else
                    matrix[u] = temp;
                u++;
            }
        }
        pushiexp(iexp);

    }

    /**
     * This function returns the Eigen vectors.
     *
     * @return the array
     */
    public double[][] getEigenVectors() {
        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }
        return Evec;
    }

    /**
     * This function returns the inverse Eigen vectors.
     *
     * @return the array
     */
    public double[][] getInverseEigenVectors() {
        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }
        return Ievc;
    }

    /**
     * This function returns the Eigen values.
     */
    public double[] getEigenValues() {
        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }
        return Eval;
    }


    // Make it a valid rate matrix (make sum of rows = 0)
    void makeValid(double[][] matrix, int dimension) {
        for (int i = 0; i < dimension; i++) {
            double sum = 0.0;
            for (int j = 0; j < dimension; j++) {
                if (i != j)
                    sum += matrix[i][j];
            }
            matrix[i][i] = -sum;
        }
    }
    
    /**
     * Complex-related methods
     * 
     */
    public void setMaxIterations(int max) {
        maxIterations = max;
    }

    public void setMaxConditionNumber(double max) {
        maxConditionNumber = max;
    }

    public void setCheckConditioning(boolean check) {
        checkConditioning = check;
    }
    
    protected void checkComplexSolutions() {
        boolean complex = false;
        for (int i = 0; i < stateCount && !complex; i++) {
            if (EvalImag[i] != 0)
                complex = true;
        }
        isComplex = complex;
    }

    public boolean getIsComplex() {
        return isComplex;
    }

    /**
     * allocate memory for the Eigen routines
     */
    protected void initialiseEigen() {

        Eval = new double[stateCount];
        Evec = new double[stateCount][stateCount];
        Ievc = new double[stateCount][stateCount];

        storedEval = new double[stateCount];
        storedEvalImag = new double[stateCount];
        storedEvec = new double[stateCount][stateCount];
        storedIevc = new double[stateCount][stateCount];
        
        amat = new double[stateCount][stateCount];

        ordr = new int[stateCount];
        evali = new double[stateCount];

        eigenInitialised = true;
        updateMatrix = true;
    }

    // Eigenvalues, eigenvectors, and inverse eigenvectors
    protected double[] Eval;
    protected double[] storedEval;
    protected double[][] Evec;
    protected double[][] storedEvec;
    protected double[][] Ievc;
    protected double[][] storedIevc;
 
    	//Imaginary Evalues
    protected double[] EvalImag;
    protected double[] storedEvalImag;
    

    List<double[][]> iexpPool = new LinkedList<double[][]>();

    private int[] ordr;
    private double[] evali;
    double amat[][];
    
    protected static final double minProb = Property.DEFAULT.tolerance();
    private static final Algebra alegbra = new Algebra(minProb);

    protected synchronized double[][] popiexp() {

        if (iexpPool.size() == 0) {
            iexpPool.add(new double[stateCount][stateCount]);
        }
        return iexpPool.remove(0);
    }

    protected synchronized void pushiexp(double[][] iexp) {
        iexpPool.add(0, iexp);
    }

	
    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

    public String toXHTML() {
        StringBuffer buffer = new StringBuffer();

        buffer.append("<em>AbsoluteCNA Model</em> (gain rate= ");
        buffer.append(gainParameter.getValue(0));
        buffer.append(", relative loss rate (*gain rate)= ");
        buffer.append(lossParameter.getValue(0));
        buffer.append(")");

        return buffer.toString();
    };
	

    
}
