package PISCA;
import java.lang.*;
import java.util.*;

import PISCA.*;
import dr.app.plugin.*;
import dr.xml.XMLObjectParser;
public class PISCALoader implements Plugin {

        public Set<XMLObjectParser> getParsers() {
                
		Set<XMLObjectParser> parsers = new HashSet<XMLObjectParser>();
//		
//		List<String> SCAparsers = Arrays.asList("AscertainedCharactersSitePatternsParser","BinarySequenceErrorModelParser","CNVSubstitutionModelParser","CenancestorTreeLikelihoodParser","LOHParser","RandomLocalClockModelParserCenancestor","RateEpochCenancestorBranchRateModelParser","RateStatisticCenancestorParser","StrictClockCenancestorBranchRatesParser"); //Add parsers here if necessary
//	
//		try {
//			for (String parserName : SCAparsers)
//			{
//								
//				parsers.add((XMLObjectParser) Class.forName(parserName).cast(Class.forName(parserName).newInstance()));
//				//parsers.add(new Class.forName("AscertainedCharactersSitePatternsParser"));
//			}
//		} catch (Exception e) {
//		}
		AscertainedCharactersSitePatternsParser parser1 = new AscertainedCharactersSitePatternsParser();
		parsers.add(parser1);
		BinarySequenceErrorModelParser parser2 = new BinarySequenceErrorModelParser();
		parsers.add(parser2);
		CNVSubstitutionModelParser parser3 = new CNVSubstitutionModelParser();
		parsers.add(parser3);
		AbsoluteCNASubstitutionModelParser parser4 = new AbsoluteCNASubstitutionModelParser();
		parsers.add(parser4);
		CenancestorTreeLikelihoodParser parser5 = new CenancestorTreeLikelihoodParser();
		parsers.add(parser5);
		LOHParser parser6 = new LOHParser();
		parsers.add(parser6);
		RandomLocalClockModelParserCenancestor parser7 = new RandomLocalClockModelParserCenancestor();
		parsers.add(parser7);
		RateEpochCenancestorBranchRateModelParser parser8 = new RateEpochCenancestorBranchRateModelParser();
		parsers.add(parser8);
		RateStatisticCenancestorParser parser9 = new RateStatisticCenancestorParser();
		parsers.add(parser9);
		StrictClockCenancestorBranchRatesParser parser10 = new StrictClockCenancestorBranchRatesParser();
		parsers.add(parser10);
		BiallelicBinarySubstitutionModelParser parser11 = new BiallelicBinarySubstitutionModelParser();
		parsers.add(parser11);
                return parsers;
        }

}
