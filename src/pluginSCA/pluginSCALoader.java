package pluginSCA;
import java.lang.*;
import java.util.*;

import pluginSCA.*;
import dr.app.plugin.*;
import dr.xml.XMLObjectParser;
public class pluginSCALoader implements Plugin {

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
		CenancestorTreeLikelihoodParser parser4 = new CenancestorTreeLikelihoodParser();
		parsers.add(parser4);
		LOHParser parser5 = new LOHParser();
		parsers.add(parser5);
		RandomLocalClockModelParserCenancestor parser6 = new RandomLocalClockModelParserCenancestor();
		parsers.add(parser6);
		RateEpochCenancestorBranchRateModelParser parser7 = new RateEpochCenancestorBranchRateModelParser();
		parsers.add(parser7);
		RateStatisticCenancestorParser parser8 = new RateStatisticCenancestorParser();
		parsers.add(parser8);
		StrictClockCenancestorBranchRatesParser parser9 = new StrictClockCenancestorBranchRatesParser();
		parsers.add(parser9);
                return parsers;
        }

}
