/*
 * BinarySequenceErrorModelParser.java
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

import dr.evolution.util.TaxonList;
import PISCA.BinarySequenceErrorModel;
import dr.inference.model.Variable;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * BinarySequenceErrorModelParser - implements a model for sequencing errors of binary data.
 *
 * @author Diego Mallo
 * @version $Id: v 1.0 $
 */
public class BinarySequenceErrorModelParser extends AbstractXMLObjectParser {

    public static final String SEQUENCE_ERROR_MODEL = "binarySequenceErrorModel";
    public static final String FN_RATE = "fnRate"; 
    public static final String FP_RATE = "fpRate";//FP == Rumen's alpha
    public static final String EXCLUDE = "exclude";
    public static final String INCLUDE = "include";
    //public static final String AGE_RELATED_RATE = "ageRelatedErrorRate";
    public static final String INDICATORS = "indicators";
 

    public String getParserName() {
        return SEQUENCE_ERROR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Variable<Double> fpRate = null;
        if (xo.hasChildNamed(FP_RATE)) {
            fpRate = (Variable) xo.getElementFirstChild(FP_RATE);
        }
        
        Variable<Double> fnRate = null;
        if (xo.hasChildNamed(FN_RATE)) {
            fnRate = (Variable) xo.getElementFirstChild(FN_RATE);
        }

        if (fpRate == null && fnRate == null) {
            throw new XMLParseException("You must specify one or other or both of " +
                    FN_RATE + " and " + FP_RATE + " parameters");
        }

        Parameter indicatorParameter = null;
        if (xo.hasChildNamed(INDICATORS)) {
            indicatorParameter = (Parameter)xo.getElementFirstChild(INDICATORS);
        }


        TaxonList includeTaxa = null;
        TaxonList excludeTaxa = null;

        if (xo.hasChildNamed(INCLUDE)) {
            includeTaxa = (TaxonList) xo.getElementFirstChild(INCLUDE);
        }

        if (xo.hasChildNamed(EXCLUDE)) {
            excludeTaxa = (TaxonList) xo.getElementFirstChild(EXCLUDE);
        }

        BinarySequenceErrorModel DamageModel = new BinarySequenceErrorModel(fpRate,fnRate,includeTaxa,excludeTaxa,indicatorParameter);

        Logger.getLogger("dr.evomodel").info("Using binary sequence error model");

        return DamageModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a model that allows for post-mortem binary-data damage.";
    }

    public Class getReturnType() {
        return BinarySequenceErrorModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
    			new ElementRule(FP_RATE, Variable.class, "The false positive likelihood", true),
            new ElementRule(FN_RATE, Variable.class, "The false negative likelihood", true),
            new ElementRule(INDICATORS, Parameter.class, "A binary indicator of whether the sequence has errors", true),
//            new XORRule(
//                    new ElementRule(INCLUDE, TaxonList.class, "A set of taxa to which to apply the damage model to"),
//                    new ElementRule(EXCLUDE, TaxonList.class, "A set of taxa to which to not apply the damage model to")
//                    , true)
    };

}
