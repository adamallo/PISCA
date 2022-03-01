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
	
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TwoStates;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
	
	
	/**
	 * LOH model
	 *
	 * @version $Id: LOH.java $
	 *
	 * @author Rumen Kostadinov
	 * @author Diego Mallo
	 */
public class LOH extends AbstractModel implements SubstitutionModel, dr.util.XHTMLable {
    
		public static final String MODEL = "model";

    	protected DataType dataType = null;

    	protected FrequencyModel freqModel;

	    private Variable<Double> alphaParameter = null;
	    private Variable<Double> betaParameter = null;

	    
	    private double alpha=0;
	    private double beta=0;
	
	    /**
	     * Constructor
	     */
	    public LOH(Variable alphaParameter, Variable betaParameter,FrequencyModel freqModel) {
	
	        //super(LOHParser.LOH_MODEL,TwoStates.INSTANCE,freqModel);
	    	super(LOHParser.LOH_MODEL);
	    	
	    	this.dataType = TwoStates.INSTANCE;

	        if (freqModel != null) {

	            if (freqModel.getDataType() != dataType) {
	                throw new IllegalArgumentException("Datatypes do not match!");
	            }

	            this.freqModel = freqModel;
	            addModel(freqModel);
	        }

	        this.alphaParameter = alphaParameter;
	        addVariable(alphaParameter);
	        alphaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
	        this.betaParameter = betaParameter;
	        addVariable(betaParameter);
	        betaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

	    }
	
	    /**
	     * set alpha
	     */
	    public void setAlpha(double alpha) {
	        alphaParameter.setValue(0, alpha);
	    }
	    
	    /**
	     * set beta
	     */
	    public void setBeta(double beta) {
	        betaParameter.setValue(0, beta);
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
	        //updateIntermediates = true; //Nothing to do
	    }

	    protected void ratesChanged() {
	        // Nothing to precalculate
	    }
	
	    /**
	     * get the complete transition probability matrix for the given distance
	     *
	     * @param distance the expected number of substitutions
	     * @param matrix   an array to store the matrix
	     */
	    public void getTransitionProbabilities(double distance, double[] matrix) {

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
	     * setup substitution matrix //Nothing to set up
	     */
	    public void setupMatrix() {
	    }
	
	    protected void setupRelativeRates() {
	    }
	
	    
	    // *****************************************************************
	    // Interface Model
	    // *****************************************************************


	    protected void storeState() {
	    } // nothing to do

	    /**
	     * Restore the additional stored state
	     */
	    protected void restoreState() {
	    } // nothing to do

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
	    }

		@Override
		public double[][] getEigenVectors() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double[][] getInverseEigenVectors() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double[] getEigenValues() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public FrequencyModel getFrequencyModel() {
			return this.freqModel;
		}

		@Override
		public DataType getDataType() {
			return this.dataType;
		}

		@Override
		protected void handleModelChangedEvent(Model model, Object object, int index) {
			// The only model registered is the frequencyModel, which does not modify this substitution model
			
		}

		@Override
		protected void handleVariableChangedEvent(Variable variable, int index, ChangeType type) {
			// We are calculating the transition probabilities on the fly, so no need to recalculate anything a priory if the rates change
			
		};
}
