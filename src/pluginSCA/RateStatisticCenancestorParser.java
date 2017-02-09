/*
 * RateStatisticParser.java
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

import pluginSCA.RateStatisticCenancestor;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.RateStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class RateStatisticCenancestorParser extends AbstractXMLObjectParser {

    public static final String RATE_STATISTIC = "rateStatisticCenancestor";
    public static final String MODE = "mode";
    public static final String MEAN = "mean";
    public static final String VARIANCE = "variance";
    public static final String COEFFICIENT_OF_VARIATION = "coefficientOfVariation";

    public String getParserName() {
        return RATE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.getAttribute(Statistic.NAME, xo.getId());
        final Tree tree = (Tree) xo.getChild(Tree.class);
        final CenancestorBranchRateModel branchRateModel = (CenancestorBranchRateModel) xo.getChild(BranchRateModel.class);

        final boolean internal = xo.getBooleanAttribute("internal");
        final boolean external = xo.getBooleanAttribute("external");

        if (!(internal || external)) {
            throw new XMLParseException("At least one of internal and external must be true!");
        }

        final String mode = xo.getStringAttribute(MODE);

        return new RateStatisticCenancestor(name, tree, branchRateModel, external, internal, mode);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns the average of the branch rates";
    }

    public Class getReturnType() {
        return RateStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(BranchRateModel.class),
            AttributeRule.newBooleanRule("internal"),
            AttributeRule.newBooleanRule("external"),
            new StringAttributeRule("mode", "This attribute determines how the rates are summarized, can be one of (mean, variance, coefficientOfVariance)", new String[]{MEAN, VARIANCE, COEFFICIENT_OF_VARIATION}, false),
            new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
    };

}
