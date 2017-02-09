/*
 * BinarySequenceErrorModel.java
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

package pluginSCA;

import dr.evolution.datatype.TwoStates;
import dr.evolution.util.TaxonList;
import pluginSCA.BinarySequenceErrorModelParser;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.model.Variable;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;

/**
 * This class incorporates uncertainty in the state at the tips of the tree and can
 * be used to model processes like sequencing error and DNA damage.
 * @author Diego Mallo
 * @version $Id$
 */
public class BinarySequenceErrorModel extends TipStatesModel {
//    public enum ErrorType {
//        TYPE_1_TRANSITIONS("type1Transitions"),
//        TYPE_2_TRANSITIONS("type2Transitions"),
//        TRANSITIONS_ONLY("transitionsOnly"),
//        ALL_SUBSTITUTIONS("allSubstitutions");
//
//
//        ErrorType(String label) {
//            this.label = label;
//        }
//
//        public String toString() {
//            return label;
//        }
//
//        final String label;
//    }
    //private final ErrorType errorType;
    private final Variable<Double> fnRate;
    private final Variable<Double> fpRate; //Rumen's alpha
    //private final Parameter ageRelatedErrorRateParameter;
    private final Parameter indicatorParameter;

    public BinarySequenceErrorModel(Variable fpRate, Variable fnRate, TaxonList includeTaxa, TaxonList excludeTaxa, Parameter indicatorParameter) {
        super(BinarySequenceErrorModelParser.SEQUENCE_ERROR_MODEL, includeTaxa, excludeTaxa);

        //this.errorType = errorType;

        if (fnRate != null) {
            this.fnRate = fnRate;
            addVariable(fnRate);
            fnRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }else{
        		this.fnRate=null;
        }
        if (fpRate !=null){
        		this.fpRate = fpRate;
        		addVariable(fpRate);
        		fpRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }else{
        		this.fpRate=null;
        }

//        if (ageRelatedErrorRateParameter != null) {
//            this.ageRelatedErrorRateParameter = ageRelatedErrorRateParameter;
//            addVariable(ageRelatedErrorRateParameter);
//        } else {
//            this.ageRelatedErrorRateParameter = null;
//        }

        if (indicatorParameter != null) {
            this.indicatorParameter = indicatorParameter;
            addVariable(indicatorParameter);
        } else {
            this.indicatorParameter = null;
        }

        if (indicatorParameter != null) {
            addStatistic(new TaxonHasErrorsStatistic());
        }
    }

    protected void taxaChanged() {
        if (indicatorParameter != null && indicatorParameter.getDimension() <= 1) {
            this.indicatorParameter.setDimension(tree.getExternalNodeCount());
        }
    }

    @Override
    public Type getModelType() {
        return Type.PARTIALS;
    }

    @Override
    public void getTipStates(int nodeIndex, int[] tipStates) {
        throw new IllegalArgumentException("This model emits only tip partials");
    }

    @Override
    public void getTipPartials(int nodeIndex, double[] partials) {

        int[] states = this.states[nodeIndex];
        if (indicatorParameter == null || indicatorParameter.getParameterValue(nodeIndex) > 0.0) {

            double pDamagedOne = 0.0;
            double pDamagedZero = 0.0;

            if (!excluded[nodeIndex]) {
                if (fnRate != null) {
                    pDamagedZero=fnRate.getValue(0);
                }
                if (fpRate !=null)
                {
                		pDamagedOne =fpRate.getValue(0);
                }

//                if (ageRelatedErrorRateParameter != null) {
//                    double rate = ageRelatedErrorRateParameter.getParameterValue(0);
//                    double age = tree.getNodeHeight(tree.getExternalNode(nodeIndex));
//                    pUndamaged *= Math.exp(-rate * age);
//                }

            }

            int k = 0;
            for (int j = 0; j < patternCount; j++) {
                switch (states[j]) {
                    case TwoStates.ZERO_STATE: // is a 0
                        partials[k] = 1-pDamagedZero; // to 0
                        partials[k + 1] = pDamagedZero; // to 1
                        break;
                    case TwoStates.ONE_STATE: // is a 1
                        partials[k] = pDamagedOne; // to 0
                        partials[k + 1] = 1-pDamagedOne; // to 1
                        break;

                    default: // is an ambiguity
                        partials[k] = 1.0;
                        partials[k + 1] = 1.0;
                }
                k += stateCount;
            }
        } else {
            int k = 0;
            for (int j = 0; j < patternCount; j++) {

                switch (states[j]) {
                    case TwoStates.ZERO_STATE: // is an A
                        partials[k] = 1.0;
                        partials[k + 1] = 0.0;

                        break;
                    case TwoStates.ONE_STATE: // is an C
                        partials[k] = 0.0;
                        partials[k + 1] = 1.0;
                        break;
                    default: // is an ambiguity
                        partials[k] = 1.0;
                        partials[k + 1] = 1.0;
                }

                k += stateCount;
            }

        }
    }
    public class TaxonHasErrorsStatistic extends Statistic.Abstract {

        public TaxonHasErrorsStatistic() {
            super("hasErrors");
        }

        public int getDimension() {
            if (indicatorParameter == null) return 0;
            return indicatorParameter.getDimension();
        }

        public String getDimensionName(int dim) {
            return taxonMap.get(dim);
        }

        public double getStatisticValue(int dim) {
            return indicatorParameter.getParameterValue(dim);
        }

    }
}
