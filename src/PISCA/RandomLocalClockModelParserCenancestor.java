/*
 * RandomLocalClockModelCenancestorParser.java
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

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.Bounds;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class RandomLocalClockModelParserCenancestor extends AbstractXMLObjectParser {

    public static final String LOCAL_BRANCH_RATES = "randomLocalClockModelCenancestor";
    public static final String RATE_INDICATORS = "rateIndicator";
    public static final String RATES = "rates";
    public static final String CLOCK_RATE = "clockRate";
    public static final String RATES_ARE_MULTIPLIERS = "ratesAreMultipliers";

    public String getParserName() {
        return LOCAL_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        Parameter rateIndicatorParameter = (Parameter) xo.getElementFirstChild(RATE_INDICATORS);
        Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES);
        Parameter meanRateParameter = null;

        if (xo.hasChildNamed(CLOCK_RATE)) {
            meanRateParameter = (Parameter) xo.getElementFirstChild(CLOCK_RATE);
        }
		
		Bounds<Double> rateBounds=meanRateParameter.getBounds();

        if (rateBounds==null|| rateBounds!=null && rateBounds.getLowerLimit(0)<0.00000001){
            Logger.getLogger("dr.evomodel").info("\n\nWARNING: If used with ascertainment bias correction, a lower-unbounded clockrate may generate numerical instability resulting in positive data likelihoods and impeding proper mixing.");
        }

        boolean ratesAreMultipliers = xo.getAttribute(RATES_ARE_MULTIPLIERS, false);

        Logger.getLogger("dr.evomodel").info("Using random local clock (RLC) model.");
        Logger.getLogger("dr.evomodel").info("  rates at change points are parameterized to be " +
                (ratesAreMultipliers ? " relative to parent rates." : "independent of parent rates."));

        return new RandomLocalClockModelCenancestor(tree, meanRateParameter, rateIndicatorParameter,
                ratesParameter, ratesAreMultipliers);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns an random local clock (RLC) model." +
                        "Each branch either has a new rate or " +
                        "inherits the rate of the branch above it depending on the indicator vector, " +
                        "which is itself sampled.";
    }

    public Class getReturnType() {
        return RandomLocalClockModelCenancestor.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new ElementRule(RATE_INDICATORS, Parameter.class, "The rate change indicators parameter", false),
            new ElementRule(RATES, Parameter.class, "The rates parameter", false),
            new ElementRule(CLOCK_RATE, Parameter.class, "The mean rate across all local clocks", true),
            AttributeRule.newBooleanRule(RATES_ARE_MULTIPLIERS, false)
    };
}
