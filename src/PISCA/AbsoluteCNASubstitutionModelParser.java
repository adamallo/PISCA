/*
 * AbsoluteCNASubstmodelParser.java
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

import dr.evolution.datatype.DataType;
import PISCA.AbsoluteCNASubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Diego Mallo
 */
public class AbsoluteCNASubstitutionModelParser extends AbstractXMLObjectParser {

	public static final String AbsoluteCNA_MODEL = "AbsoluteCNAModel";
	public static final String GAIN_RATE = "gain_rate";
    public static final String RELATIVE_LOSS_RATE= "relative_loss_rate";
    public static final String FREQUENCIES = "frequencies";
	public static final String MAX_STATE = "maxState";

    public String getParserName() {
        return AbsoluteCNA_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Variable gainParam = (Variable) xo.getElementFirstChild(GAIN_RATE);
        Variable lossParam = (Variable) xo.getElementFirstChild(RELATIVE_LOSS_RATE);
        //FrequencyModel freqModel = (FrequencyModel) xo.getElementFirstChild(FrequencyModelParser.FREQUENCIES);
        XMLObject cxo = xo.getChild(FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);
        DataType dataType = freqModel.getDataType();
        
        Logger.getLogger("dr.evomodel").info("Creating AbsoluteCNA substitution model. Initial gain_rate = " +
                gainParam.getValue(0) + ", relative_loss_rate = " + lossParam.getValue(0));

        return new AbsoluteCNASubstitutionModel(dataType,gainParam,lossParam,freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an instance of the AbsoluteCNA (Mallo, Malley 2022) model of SGA evolution.";
    }

    public Class getReturnType() {
        return AbsoluteCNASubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
    		new ElementRule(FREQUENCIES, FrequencyModel.class),
            new ElementRule(GAIN_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}),
            new ElementRule(RELATIVE_LOSS_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}),
    };
}
