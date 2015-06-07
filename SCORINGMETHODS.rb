require "DNATOOLBOX.rb"
require "HASHFILE.rb";
require "HOME.rb";
require "GZ.rb";
NaN="NaN"; #NaN="NaN"; #changed 20130625
class Float
	def to_s()
		return sprintf("%g",self);
	end
end
class ScoringFtn
	def initialize()
	end
	def length()
		raise ("UNIMPLEMENTED FTN length()" );
	end
	def score(seq,pos) #TODO: zero indexed?
		raise ("UNIMPLEMENTED FTN score(seq,pos)" );
	end
	def rcScore(seq,pos)
		raise ("UNIMPLEMENTED FTN rcScore(seq,pos)" );
	end
end

class TestScoringFtn <ScoringFtn
	def initialize()
		@temp=0;
	end
	def length()
		return 1;
	end
	def score(seq, pos)
		@temp+=1;
		return @temp;
	end
	def rcScore(seq, pos)
		@temp-=1;
		return @temp;
	end

end
	
class LookupScoringFtn < ScoringFtn
	def initialize(filePath)
		@scores=Hash.new();
		lastLen=-1;
		File.foreach(filePath){|line|
			line.chomp!();
			next if line.nil? || line=="" || line[0,1]=="#";
			
			seq, theScore = line.split("\t");
			theScore = theScore.to_f();
			if lastLen==-1
				lastLen=seq.length();
			elsif lastLen!=seq.length()
				raise("Multiple lengths!!! in "+ filePath+"\n");
			end
			@scores[seq]=theScore;
		}
		@myLength=lastLen;
		@myFilePath = filePath;
	end
	def length()
		return @myLength;
	end
	def score(seq, pos)
		if self.length!=seq[pos,self.length].length
			raise("ERROR!  Checking too far!!! "+pos.to_s+"\n");
		end
		if seq[pos,self.length]=~/[MRWSYKVHDBXN]/
			theScore = 0.0;
			count = 0;
			#expand
			expandIUPAC(seq[pos,self.length]).each{|posSeq|
				if @scores[posSeq].nil?
					if @scores[revcomp(posSeq)].nil?
						raise("CAN'T Find poss seq "+posSeq+" from "+seq[pos,self.length]+" in hash "+@myFilePath+"\n");
					end
					theScore+= @scores[revcomp(posSeq)];
				else
					theScore+=@scores[posSeq];
				end
				count+=1;
			}
			return theScore/count;
		end
		if @scores[seq[pos,self.length]].nil?
			if @scores[revcomp(seq[pos,self.length])].nil?
				raise("CAN'T Find "+seq[pos,self.length]+" in hash "+@myFilePath+"\n");
			end
			return @scores[revcomp(seq[pos,self.length])];
		end
		return @scores[seq[pos,self.length]];
	end
	def rcScore(seq, pos)
		if self.length!=seq[pos,self.length].length
			raise("ERROR!  Checking too far!!! "+pos.to_s+"\n");
		end
		if seq[pos,self.length]=~/[MRWSYKVHDBXN]/
			theScore = 0.0;
			count = 0;
			#expand
			expandIUPAC(seq[pos,self.length]).each{|posSeq|
				if @scores[revcomp(posSeq)].nil?
					if @scores[posSeq].nil?
						raise("CAN'T Find poss seq "+posSeq+" from "+seq[pos,self.length]+" in hash "+@myFilePath+"\n");
					end
					theScore+= @scores[posSeq];
				else
					theScore+=@scores[revcomp(posSeq)];
				end
				count+=1;
			}
			return theScore/count;
		end
		if @scores[revcomp(seq[pos,self.length])].nil?
			if @scores[seq[pos,self.length]].nil?
				raise("CAN'T Find "+seq[pos,self.length]+" in hash "+@myFilePath+"\n");
			end
			return @scores[seq[pos,self.length]];
		end
		return @scores[revcomp(seq[pos,self.length])];
	end
end

class ReverseLookupScoringFtn < LookupScoringFtn
	alias_method :score, :rcScore;
	alias_method :rcScore, :score;
	#def rcScore(seq, pos)
	#	LookupScoringFtn.instance_method(:score).bind(self).call(seq,pos);
	#end
	#def score(seq,pos)
	#end
end
class ConsensusScoringFtn < ScoringFtn
	def initialize(motif)
		allPossibilities = expandIUPAC(motif);
		@isOne = Hash.new();
		@rcIsOne=Hash.new();
		allPossibilities.each{|possib|
			@isOne[possib]=1;
			@rcIsOne[revcomp(possib)]=1;
		}
		@myLength=motif.length();
	end
	def length()
		return @myLength;
	end
	def score(seq, pos)
		if self.length!=seq[pos,self.length].length
			raise("ERROR!  Checking too far!!! "+pos.to_s+" "+seq+"\n");
		end
		if @isOne[seq[pos,self.length]]!=nil
			return 1;
		end
		return 0;
	end
	def rcScore(seq, pos)
		if @rcIsOne[seq[pos,self.length]]!=nil
			return 1;
		end
		return 0;
	end
end

class GCContentScoringFtn < ConsensusScoringFtn
	def initialize()
		super("S");
	end
end

class AContentScoringFtn < ConsensusScoringFtn
	def initialize()
		super("A");
	end
end
class TContentScoringFtn < ConsensusScoringFtn
	def initialize()
		super("T");
	end
end
class GContentScoringFtn < ConsensusScoringFtn
	def initialize()
		super("G");
	end
end
class CContentScoringFtn < ConsensusScoringFtn
	def initialize()
		super("C");
	end
end

class NonSSConsensusScoringFtn < ConsensusScoringFtn
	def initialize(motif)
		super(motif)
		@rcIsOne.keys.each{|key|
			@isOne[key]=1;
		}
		@rcIsOne=@isOne;
	end
end
class NonSSPWMScoringFtn <ScoringFtn
	def initialize()
	end
	def length()
		return @PWM['A'].length();
	end
	
	def score(seq, pos)
		curScore = 0.0;
		0.upto(self.length-1){|j|
			curScore+=@PWM[seq[pos+j,1]][j];
		}
		rScore = 0.0;
		0.upto(self.length-1){|j|
			rScore+=@rcPWM[seq[pos+j,1]][j];
		}
		return [curScore, rScore].max();
	end
	def rcScore(seq,pos)
		return score(seq,pos);
	end

end

class NonSSPLPMFromFileScoringFtn < NonSSPWMScoringFtn
	def initialize(filePath)
		thePFM = smartLoadPFMFromFile(filePath);
		@PWM = thePFM.toPLPM();
		@rcPWM = thePFM.revcomp().toPLPM();
	end
end

class NonSSPWMFromFileScoringFtn < NonSSPWMScoringFtn
	def initialize(filePath)
		thePFM = smartLoadPFMFromFile(filePath);
		@PWM = thePFM.toPWM();
		@rcPWM = thePFM.revcomp().toPWM();
	end
end

class NonSSPWMFromIUPACScoringFtn < NonSSPWMScoringFtn
	def initialize(iupacStr)
		iupacHash = loadHashFile("IUPAC_standard_and_lowC_no0s.txt");
		thePFM = PFM.new(nil);
		thePFM.loadFromIUPAC(iupacStr, iupacHash);
		@PWM = thePFM.toPWM();
		@rcPWM = thePFM.revcomp().toPWM();
	end
end

class PWMScoringFtn < ScoringFtn
	def initialize()
	end
	def length()
		return @PWM['A'].length();
	end
	
	def score(seq, pos)
		curScore = 0.0;
		0.upto(self.length-1){|j|
			if @PWM[seq[pos+j,1]]!=nil
				curScore+=@PWM[seq[pos+j,1]][j];
			end
		}
		return curScore;
	end
	def rcScore(seq,pos)
		curScore = 0.0;
		  0.upto(self.length-1){|j|
			if @rcPWM[seq[pos+j,1]]!=nil
				curScore+=@rcPWM[seq[pos+j,1]][j];
			end
		  }
		return curScore;
	end
end


class NeutralPWMFromFileScoringFtn < PWMScoringFtn
	def initialize(filePath)
		thePFM = smartLoadPFMFromFile(filePath);
		@PWM = thePFM.toPWM({"A"=>0.25, "T"=>0.25, "G"=>0.25, "C"=>0.25});
		@rcPWM = thePFM.revcomp().toPWM({"A"=>0.25, "T"=>0.25, "G"=>0.25, "C"=>0.25});
	end
end

class PLPMFromFileScoringFtn < PWMScoringFtn
	def initialize(filePath)
		thePFM = smartLoadPFMFromFile(filePath);
		@PWM = thePFM.toPLPM();
		@rcPWM = thePFM.revcomp().toPLPM();
	end
end

class PWMFromFileScoringFtn < PWMScoringFtn
	def initialize(filePath)
		thePFM = smartLoadPFMFromFile(filePath);
		@PWM = thePFM.toPWM();
		@rcPWM = thePFM.revcomp().toPWM();
	end
end
class ONLYBESTPWMFromFileScoringFtn < PWMFromFileScoringFtn
  def score(seq,pos)
		score = super(seq,pos);
		if score>0
			return score
		else
			return 0;
		end
	end
  def rcScore(seq,pos)
		score = super(seq,pos);
		if score>0
			return score
		else
			return 0;
		end
	end
end
class PWMFromIUPACScoringFtn < PWMScoringFtn
	def initialize(iupacStr)
		iupacHash = loadHashFile("IUPAC_standard_and_lowC_no0s.txt");
		thePFM = PFM.new(nil);
		thePFM.loadFromIUPAC(iupacStr, iupacHash);
		@PWM = thePFM.toPWM();
		@rcPWM = thePFM.revcomp().toPWM();
	end
end

#	Sequence Scorer is used to score a region of DNA, not a single position
class SequenceScorer
	def initialize(myScoringFunction)
		@scoreFtn = myScoringFunction;
	end
	def length()
		return @scoreFtn.length();
	end
	def scoreThisSeq(seq, stPos, enPos, strand)
		raise ("UNIMPLEMENTED FTN scoreThisSeq()" );
	end
end

class ExternalPreCalcScorer <SequenceScorer
	def scoreThisSeq(seq, stPos, enPos, strand, chr)
		raise ("UNIMPLEMENTED FTN scoreThisSeq()" );
	end
end

class WigScorer<ExternalPreCalcScorer
	require "GENOMEDATA.rb";
	def initialize(wigFPs)
		fs = wigFPs.split(":");
		if fs.length()==1
			@myWigPlus = readWigFile(nil,fs[0]);
			@myWigMinus = @myWigPlus;
		elsif fs.length==2
			@myWigPlus = readWigFile(nil,fs[0]);
			@myWigMinus = readWigFile(nil,fs[1]);
		else
			raise("Wrong number of args (#{wigFPs}) for WigScorer");
		end
	end
	
	def howToCombine(scores)
		raise("UNIMPLEMENTED FTN howToAdd(oldScore, curScore, curPUnP)" );
	end
	def scoreThisSeq(seq, stPos, enPos, strand, chr)
		if strand=="+"
			return howToCombine(@myWigPlus[chr][stPos,enPos-stPos+1]);
		else
			return howToCombine(@myWigMinus[chr][stPos,enPos-stPos+1]);
		end
	end
end

class SumWigScorer < WigScorer
	def howToCombine(scores)
		return scores.inject(0.0) { |sum, el| sum + el };
	end
end
class AverageWigScorer < WigScorer
	def howToCombine(scores)
		return scores.inject(0.0) { |sum, el| sum + el } / scores.length();
	end
end


class ExternalPreCalc
	def initialize(pUnMapFP)
		#load the unpaired map
		@accFiles =Hash.new();
		@pUnP = Hash.new();
		File.foreach(pUnMapFP){|line|
			chr, path = line.chomp.split("\t");
			@accFiles[chr+"+"] = path+".f.plout";
			@accFiles[chr+"-"] = path+".r.plout";
		}
	end
	def loadFile(file)
		data =[];
		File.foreach(file){|sline|
			sline.chomp!();
			if sline==nil||sline==""||sline[0,1]=="#"
				next;
			end
			pos, probUnpaired = sline.split(" ");
			data.push(probUnpaired.to_f());
		}
		return data;
	end
	def [](key)
		if @pUnP[key].nil?
			@pUnP[key]= loadFile(@accFiles[key]);
		end
		return @pUnP[key];
	end
end


class  PUnpairedScorer< ExternalPreCalcScorer
	def initialize(myScoringFunction)
		raise ("Need 2 parameters for this scorer." );
	end
	def initialize(myScoringFunction, pUnMap)
		super(myScoringFunction);
		@preCalc = pUnMap;
	end
	
	def howToAdd(oldScore, curScore, curPUnP)
		raise ("UNIMPLEMENTED FTN howToAdd(oldScore, curScore, curPUnP)" );
	end
	def finalTreatment(totalScore, count)
		raise ("UNIMPLEMENTED FTN finalTreatment(totalScore, count)" );
	end
	def scoreThisSeq(seq, stPos, enPos, strand, chr)
		totalScore=0.0;
		enPos = enPos-(@scoreFtn.length-1);
		count=0;
		stPos.upto(enPos){|curStart|
			if strand=="+"
				curScore = @scoreFtn.score(seq, curStart);
			else
				curScore= @scoreFtn.rcScore(seq, curStart);
			end
			thisPUnP = 1.0;
			curStart.upto(curStart+@scoreFtn.length()-1){|curPos|
				if @preCalc[chr+strand][curPos].nil?
					print("Warning: going past chromosome end for PunP "+[stPos, enPos, strand, chr, curPos].join(" ")+"\n")
				else
					thisPUnP = @preCalc[chr+strand][curPos];
				end
			}
			totalScore = howToAdd(totalScore, curScore, thisPUnP);
			count+=1;
		}
		return finalTreatment(totalScore, count);
	end
end

class AdditivePUnpairedScorer < PUnpairedScorer
	def howToAdd(oldTotal, curScore, curPUnP)
		if @scoreFtn.is_a?(PWMScoringFtn)
			totalScore = curScore+ Math.log(thisPUnP)/Math.log(2);
			if totalScore>0
				newTotal = oldTotal+totalScore;
			else
				newTotal= oldTotal + curScore+ Math.log(thisPUnP)/Math.log(2);
			end
		else
			if curPUnP.nil?
				raise("bad PUnP: "+[oldTotal, curScore, curPUnP].join(" ")+"\n");
			elsif curScore.nil?
				raise("bad curScore: "+[oldTotal, curScore, curPUnP].join(" ")+"\n");
			end
			newTotal= oldTotal+ (curScore*curPUnP);
		end
		return newTotal;
	end
	def finalTreatment(totalScore, count)
		return totalScore/count;
	end

end

class BestPUnpairedScorer < PUnpairedScorer
	def howToAdd(oldTotal, curScore, curPUnP)
		if @scoreFtn.is_a?(PWMScoringFtn)
			newScore= curScore+ Math.log(thisPUnP)/Math.log(2);
		else
			newScore =  curScore*curPUnP;
		end
		return [oldTotal, newScore].max();
	end
	def finalTreatment(totalScore, count)
		return totalScore;
	end
end

class HairpinLoopScorer < SequenceScorer
	def initialize(myScoringFtn)
		super(myScoringFtn);
		@fHash=Hash.new();
		@bHash=Hash.new();
		@rnaBPFEs = [];
		firstLine = true;
		lineNum=0;
		File.foreach(HOME+"/lisensing/RNA_dinucFEs.txt"){|line|
			if firstLine
				firstLine=false;
				theHeader = line.chomp.split("\t");
				theHeader.shift();
				0.upto(theHeader.length-1){|i|
					@bHash[theHeader[i]]=i;
				}
				next;
			end

			theData = line.chomp.split("\t");
			@fHash[theData.shift()]=lineNum;
			theData.map!{|v|v.to_f};
			@rnaBPFEs.push(theData);
			lineNum+=1;
		}
		@getThisMuch=8;
	end
	def getFE(seq1, seq2)
		if seq1=~/[MRWSYKVHDBNX]/ || seq2=~/[MRWSYKVHDBNX]/
			sumFE = 0.0;
			count=0;
			s1E = expandIUPAC(seq1);
			s2E = expandIUPAC(seq2);
			s1E.each{|s1ex|
				s2E.each{|s2ex|
					sumFE+=@rnaBPFEs[@fHash[s1ex]][@bHash[s2ex]];
					count+=1;
				}
			}
			return sumFE/count;
		end
		return(@rnaBPFEs[@fHash[seq1]][@bHash[seq2]]);
	end
	def howToAdd(thisNFE, total)
		return [thisNFE,total].max();
	end
	def finalTreatment(totalScore,count)
		return(totalScore/count);
	end
	def scoreThisSeq(seq, stPos, enPos, strand)
		totalScore=0.0;
		enPos = enPos-(@scoreFtn.length-1);
		count=0;
		stPos.upto(enPos){|curStart|
			if strand=="+"
				curScore = @scoreFtn.score(seq, curStart);
			else
				curScore= @scoreFtn.rcScore(seq, curStart);
			end
			if curScore>0 && (curStart-@getThisMuch)>0 && (curStart+@scoreFtn.length+@getThisMuch)<seq.length()#hit
				#take 8nt to either side and see if they line up
				front = seq[curStart+@scoreFtn.length, @getThisMuch];
				back = seq[curStart-@getThisMuch, @getThisMuch].reverse();
				totFE = 0;
				leastFE = 0;
				0.upto(@getThisMuch-2){|i|
					curFE = self.getFE(front[i,2], back[i,2]);
					if curFE>0 && i==0
						break;
					end
					totFE+=curFE;
					leastFE = [totFE, leastFE].min();
				}
				totalScore=howToAdd(-leastFE, totalScore);
			end
			count+=1;
		}
		return finalTreatment(totalScore, count);
	end

end

class RatioScorer < SequenceScorer
	def initialize(scorer1, scorer2)
		@scorer1=scorer1;
		@scorer2=scorer2;
	end
	def length()
		return [@scorer1.length(), @scorer2.length()].max();
	end
	def scoreThisSeq(seq, stPos, enPos, strand) #coordinates in the string, start from 0, inclusive
		if @scorer1.scoreThisSeq(seq, stPos, enPos, strand)==0 && @scorer2.scoreThisSeq(seq, stPos, enPos, strand)==0
			ratio = 1;	
		elsif @scorer1.scoreThisSeq(seq, stPos, enPos, strand)==0
			ratio = 0.5/(enPos-stPos+1);
		elsif @scorer2.scoreThisSeq(seq, stPos, enPos, strand)==0
			ratio = (enPos-stPos+1)/0.5;
		else
			ratio = @scorer1.scoreThisSeq(seq, stPos, enPos, strand)/@scorer2.scoreThisSeq(seq, stPos, enPos, strand);
		end
		return Math.log(ratio);
	end
end

class PolyNScorer < SequenceScorer #Fondufe-Mittendorf, 2008
	def initialize(n)
		@myN=n;
		@myRCN = revcomp(n);
	end
	def scoreThisSeq(seq, stPos, enPos, strand) #coordinates in the string, start from 0, inclusive
		totalScore=0.0;
		count=0;
		roleTop = 0;
		roleBottom=0;
		stPos.upto(enPos){|curStart|
			if seq[curStart,1]==@myN
				roleTop+=1;
				if roleTop>=5
					if roleTop==5
						totalScore+=1;
					else
						totalScore+=0.2;
					end
				end
				roleBottom=0;
			elsif seq[curStart,1]==@myRCN
				roleBottom+=1;
				if roleBottom>=5
					if roleBottom==5
						totalScore+=1;
					else
						totalScore+=0.2;
					end
				end
				roleTop=0;
			else 
				roleBottom=0;
				roleTop=0;
			end
			count+=1;
		}
		return totalScore/count;
	end
	def length()
		return 1;
	end
end

class SSPolyNScorer < PolyNScorer #Fondufe-Mittendorf, 2008
	def scoreThisSeq(seq, stPos, enPos, strand) #coordinates in the string, start from 0, inclusive
		totalScore=0.0;
		count=0;
		roleTop = 0;
		roleBottom=0;
		stPos.upto(enPos){|curStart|
			if strand=="+" && seq[curStart,1]==@myN #scanning top strand
				roleTop+=1;
				if roleTop>=5
					if roleTop==5
						totalScore+=1;
					else
						totalScore+=0.2;
					end
				end
				roleBottom=0;
			elsif strand=="-" && seq[curStart,1]==@myRCN
				roleBottom+=1;
				if roleBottom>=5
					if roleBottom==5
						totalScore+=1;
					else
						totalScore+=0.2;
					end
				end
				roleTop=0;
			else 
				roleBottom=0;
				roleTop=0;
			end
			count+=1;
		}
		return totalScore/count;
	end
end

class GQuad < SequenceScorer #Capra, 2010
	def initialize()
	end
	def scoreThisSeq(seq, stPos, enPos, strand) #coordinates in the string, start from 0, inclusive
		if seq[stPos, enPos-stPos+1]=~/GGG[ATGC]{2,25}GGG[ATGC]{2,25}GGG[ATGC]{2,25}GGG/
			return 1;
		elsif seq[stPos, enPos-stPos+1]=~/CCC[ATGC]{2,25}CCC[ATGC]{2,25}CCC[ATGC]{2,25}CCC/#on bottom
			return 1;
		else
			return 0;
		end
	end
	def length()
		return 1;
	end
end

class SimpleScorer < SequenceScorer
	def finalOperation(score, count)
		raise("NOT IMP: finalOperation(s,c) SimpleScorer");
	end
	def scoreThisSeq(seq, stPos, enPos, strand) #coordinates in the string, start from 0, inclusive
		totalScore=0.0;
		enPos = enPos-(@scoreFtn.length-1);
		count=0;
		if strand=="+"
			stPos.upto(enPos){|curStart|
				totalScore+= @scoreFtn.score(seq, curStart);
				count+=1;
			}
		else
			stPos.upto(enPos){|curStart|
				totalScore+= @scoreFtn.rcScore(seq, curStart);
				count+=1;
			}
		end
		return finalOperation(totalScore,count);
	end
end

class SimpleScorerSym < SimpleScorer
	def scoreThisSeq(seq, stPos, enPos, strand) #coordinates in the string, start from 0, inclusive
		totalScore=0.0;
		enPos = enPos-(@scoreFtn.length-1);
		count=0;
		#include both top and bottom strands
		stPos.upto(enPos){|curStart|
			totalScore+= @scoreFtn.score(seq, curStart);
			count+=1;
		}
		stPos.upto(enPos){|curStart|
			totalScore+= @scoreFtn.rcScore(seq, curStart);
			count+=1;
		}
		return finalOperation(totalScore,count);
	end
end
class SimpleSum <SimpleScorer
	def finalOperation(score, count)
		return (score);
	end
end
class SimpleAverage <SimpleScorer
	def finalOperation(score, count)
		return (score/count);
	end
end

class SimpleSumSym <SimpleScorerSym
	def finalOperation(score, count)
		return (score);
	end
end
class SimpleAverageSym <SimpleScorerSym
	def finalOperation(score, count)
		return (score/count);
	end
end

class MultByConstant <SequenceScorer
	def initialize( baseScorer, constant)
		@myConst = constant;
		@myBaseScorer = baseScorer;
	end
	def length()
		return @myBaseScorer.length()
	end
	def scoreThisSeq(seq, stPos, enPos, strand)
		return @myBaseScorer.scoreThisSeq(seq, stPos, enPos, strand)*@myConst;
	end
end

class BendScorer <MultByConstant
	def initialize(scoringFtn)
		super(SimpleSum.new(scoringFtn), (3.46543972*(10**-17)));#(1.0*(10**-9)*3.46543972*(10**-8)));
	end
end

class MeltingTemp < SequenceScorer
	def initialize()
		@scoreHash = Hash.new();
		@scoreHash["AA"]=[-7.9, -22.2];
		@scoreHash["TT"]=[-7.9, -22.2];
		@scoreHash["AT"]=[-7.2, -20.4];
		@scoreHash["TA"]=[-7.2, -21.3];
		@scoreHash["CA"]=[-8.5, -22.7];
		@scoreHash["AC"]=[-8.4, -22.4];
		@scoreHash["GT"]=[-8.4, -22.4];
		@scoreHash["TG"]=[-8.5, -22.7];
		@scoreHash["CT"]=[-7.8, -21];
		@scoreHash["TC"]=[-8.2, -22.2];
		@scoreHash["GA"]=[-8.2, -22.2];
		@scoreHash["AG"]=[-7.8, -21];
		@scoreHash["CG"]=[-10.6, -27.2];
		@scoreHash["GC"]=[-9.8, -24.4];
		@scoreHash["GG"]=[-8, -19.9];
		@scoreHash["CC"]=[-8, -19.9];
	end
	def length()
		return 2;
	end
	def scoreThisSeq(seq, stPos, enPos, strand)
		enPos = enPos-(2-1);
		count=0;
		hTot=0.0;
		sTot=0.0;
		hitGCCG=false;
		stPos.upto(enPos){|curStart|
			if seq[curStart,2]=~/[MRWSYKVHDBNX]/
				sumH = 0.0;
				sumS = 0.0;
				countExpand = 0;
				expandIUPAC(seq[curStart,2]).each{|posSeq|
					curH, curS = @scoreHash[posSeq];
					sumH+=curH;
					sumS+=curS;
					countExpand+=1;
				}
				h=sumH/countExpand;
				s=sumS/countExpand;
			else
				h,s = @scoreHash[seq[curStart,2]];
			end
			hTot+=h;
			sTot+=s;
			count+=1;
			if seq[curStart,2]=="GC" || seq[curStart,2]=="CG"
				hitGCCG=true;
			end
		}
		if hitGCCG
			hTot+=0.1;
			sTot-=2.8;
		else
			hTot+=2.3;
			sTot+=4.1;
		end
		hTot = hTot*1000;
		
		return (hTot/(sTot-32.96)) - 273.5;
	end
end

class ThresholdAdditive < SequenceScorer
	def scoreThisSeq(seq, stPos, enPos, strand) #coordinates in the string, start from 0, inclusive
		totalScore=0.0;
		enPos = enPos-(@scoreFtn.length-1);
		isTop = strand=="+";
		stPos.upto(enPos){|curStart|
			if isTop
				curScore = @scoreFtn.score(seq, curStart);
			else
				curScore = @scoreFtn.rcScore(seq, curStart);
			end
			if curScore>0
				totalScore+=curScore;
			end
		}
		return totalScore;
	end
end

class SequenceProfiler
	def initialize(myScoringFtn)
		@myScoringFtn = myScoringFtn;
	end
	def makeProfile(seq) 
		theProfile = [];
		0.upto(seq.length-@myScoringFtn.length()){|curStart|
			theProfile.push(@myScoringFtn.score(seq, curStart));
		}
		return theProfile;	
	end
end

class GFFProfiler
	def initialize(gffFP, genomeFP, myScoringFtn, flankingSeq)
		@myScoringFtn = myScoringFtn;
		@allGFFs = [];
		count=0;
		totalLength=0;
		@maxFeatureLength=0;
		File.foreach(gffFP){|line|
			line.chomp!();
			if line==nil||line==""||line[0,1]=="#"||line[0,1]==">"||line[0,11]=="track name="
				next;
			end
			allData = line.split("\t");
			#chr start end strand name
			st = allData[3].to_i()-1;
			en = allData[4].to_i()-1;
			if st>en
				warn("Start is greater than end for GFF entry, so I will swap them: #{line}\n");
				temp=st;
				st=en;
				en=temp;
			end
			fLen=en-st+1;
			@maxFeatureLength=[@maxFeatureLength,fLen].max();
			@allGFFs.push([allData[0], st, en, allData[6], allData[8]]);
			totalLength+=  fLen;#-@myScoringFtn.length+1;
			count+=1;
		}
		@avgFeatureLength = totalLength/count; 
		#print(count.to_s()+" "+ totalLength.to_s()+" "+@avgFeatureLength.to_s()+"\n");
		@chromosomes =ChrSeq.new(genomeFP);
		@flank=flankingSeq;
		@useAvg=false;
		@scaleCentres=true;
	end
	
	def revGFFs()
		0.upto(@allGFFs.length()-1){|i|
			if @allGFFs[i][3]=="+"
				@allGFFs[i][3]="-";
			else
				@allGFFs[i][3]="+";
			end
		}
	end
	def makeNaNProfile()
		return [NaN]*(@avgFeatureLength+@flank*2);
	end

	def useAvg()
		@useAvg=true;
	end
	def dontScale();
		@scaleCentres=false;
	end

	def reshapeData(curProfile,curLen, cPIndex)
		if @useAvg
			return reshapeDataAvg(curProfile,curLen, cPIndex);
		else
			return reshapeDataSum(curProfile,curLen, cPIndex);
		end
	end

	def reshapeDataAvg(curProfile,curLen, cPIndex)
		#combines data by average
		thisMid = [0.0]*@avgFeatureLength;
		theCounts = [0.0]*@avgFeatureLength;
		mapRatio=(0.0+@avgFeatureLength)/curLen;
		if curLen>=@avgFeatureLength#compress
			0.upto(curLen-1){|j|
				dest=(mapRatio*(j));
				ceilPct =0;
				floor=dest.floor();
				ceil=(dest+mapRatio).floor();
				if ceil>floor
					ceilPct=(dest+mapRatio- (dest+mapRatio).floor())/mapRatio;
				end
				floorPct=1-ceilPct;
				#print([curLen, @avgFeatureLength, dest, mapRatio, floorPct, ceilPct, floor, ceil].join("\t")+"\n");
				#p [cPIndex, j, floorPct, floor];
				#p curProfile;
				next if curProfile[cPIndex+j]==NaN || curProfile[cPIndex+j].nil?; #skip if data has a NaN

				thisMid[floor]+=floorPct*curProfile[cPIndex+j];
				theCounts[floor]+=floorPct;
				if ceil<@avgFeatureLength
					thisMid[ceil]+=ceilPct*curProfile[cPIndex+j];
					theCounts[ceil]+=ceilPct;
				end
			}
		else#expand
			0.upto(curLen-1){|j|
				dest = (mapRatio*(j));
				nextDest = dest+mapRatio;
				first=dest.floor();
				last=(nextDest).floor()
				totalPct=0.0;
				first.upto(last){|l|
					if l==@avgFeatureLength
						next;
					elsif l==first
						curPct = (1-(dest-dest.floor()))/mapRatio;
					elsif l==last
						curPct = (nextDest-nextDest.floor())/mapRatio;
					else
						curPct=1.0/mapRatio;
					end
					#print([curLen, @avgFeatureLength, dest, nextDest, mapRatio, first, last, l, curPct].join("\t")+"\n");
					begin
						if !(curProfile[cPIndex+j].nil? || curProfile[cPIndex+j]==NaN)
							theCounts[l]+=curPct;
							thisMid[l]+=curPct*curProfile[cPIndex+j];
						end
					rescue TypeError
						p(thisMid[l]);
						p(curPct);
						p(cPIndex);
						p(curProfile[cPIndex+j]);
						raise();
					end
					totalPct+=curPct;
				}
				if totalPct-1.0>=0.001
					print("ERROR: Total Percent is not 1: "+totalPct.to_s()+"\n");
				end
			}
		end
		0.upto(thisMid.length()-1){|i|
			if theCounts[i]==0
				thisMid[i]=NaN;
			else
				thisMid[i] = thisMid[i]/theCounts[i];
			end
		}
		return thisMid;
	end

	def reshapeDataSum(curProfile,curLen, cPIndex)
		#combines data by sum -unsure if this is actually working correctly because it's a weird method - I think it mght also average
		thisMid = [0.0]*@avgFeatureLength;
		mapRatio=(0.0+@avgFeatureLength)/curLen;
		if curLen>=@avgFeatureLength#compress
			0.upto(curLen-1){|j|
				dest=(mapRatio*(j));
				ceilPct =0;
				floor=dest.floor();
				ceil=(dest+mapRatio).floor();
				if ceil>floor
					ceilPct=(dest+mapRatio- (dest+mapRatio).floor())/mapRatio;
				end
				floorPct=1-ceilPct;
				#print([curLen, @avgFeatureLength, dest, mapRatio, floorPct, ceilPct, floor, ceil].join("\t")+"\n");
				#p [cPIndex, j, floorPct, floor];
				#p curProfile;
				next if curProfile[cPIndex+j]==NaN; #skip if data has a NaN
				#p([curProfile.length(),j,cPIndex, curLen]);
				thisMid[floor]+=floorPct*curProfile[cPIndex+j];
				if ceil<@avgFeatureLength
					thisMid[ceil]+=ceilPct*curProfile[cPIndex+j];
				end
			}
		else#expand
			0.upto(curLen-1){|j|
				dest = (mapRatio*(j));
				nextDest = dest+mapRatio;
				first=dest.floor();
				last=(nextDest).floor()
				totalPct=0.0;
				first.upto(last){|l|
					if l==@avgFeatureLength
						next;
					elsif l==first
						curPct = (1-(dest-dest.floor()))/mapRatio;
					elsif l==last
						curPct = (nextDest-nextDest.floor())/mapRatio;
					else
						curPct=1.0/mapRatio;
					end
					#print([curLen, @avgFeatureLength, dest, nextDest, mapRatio, first, last, l, curPct].join("\t")+"\n");
					begin
					if curProfile[cPIndex+j]!=NaN
						thisMid[l]+=curPct*curProfile[cPIndex+j];
					end
					rescue TypeError
						p(thisMid[l]);
						p(curPct);
						p(cPIndex);
						p(curProfile[cPIndex+j]);
						raise();
					end
					totalPct+=curPct;
				}
				if totalPct-1.0>=0.001
					print("ERROR: Total Percent is not 1: "+totalPct.to_s()+"\n");
				end
			}
		end
		return thisMid;
	end
	
	def doThisGFF(chr, st, en, strand)
		if st> @chromosomes[chr].length()-1|| en <0
			raise("Entry beyond end of chromosome "+chr+" "+st.to_s+" "+en.to_s+" "+strand+"\n");
		end
		if strand=="-"
			strandOffsetSt = @myScoringFtn.length()-1;
			strandOffsetEnd=0;
		else
			strandOffsetEnd = @myScoringFtn.length()-1;
			strandOffsetSt=0;
		end
		#print([chr, st, en, strand].join(" "));
		idealStart=st-@flank-strandOffsetSt;
		fStart = [0,idealStart].max();
		idealEnd=en+@flank-strandOffsetSt;
		fEnd = [@chromosomes[chr].length()-1, idealEnd].min();
		#fEnd-= (@myScoringFtn.length()-1);
		#print(" "+[fStart, fEnd, @chromosomes[chr].length].join(" ")+"\n");
		if st<0 || en>@chromosomes[chr].length()-1
			print("skipping this entry due to going off end of chr.\n");
			return makeNaNProfile();
		end
		curProfile = self.makeProfile(chr, fStart, fEnd, strand);
		fullProfile = [];
		#now add NaNs to front if appropriate

		#add initial NaNs
		(idealStart).upto(fStart-1){|i| #NaNs at start
			fullProfile.push(NaN);
		}
		cPIndex=0;
		
		if @scaleCentres
			(fStart).upto(st-1-strandOffsetSt){|i| #non-NaN start
				fullProfile.push(curProfile[cPIndex]);
				cPIndex+=1;
			}
			#now scale middle bin accordingly to @avgFeatureLength.
			curLen=(en - st+1); 
			#p([curProfile.length(),curLen,cPIndex]); #exit;
			thisMid = reshapeData(curProfile, curLen, cPIndex);
			fullProfile.concat(thisMid.map{|e|e.to_s});
			cPIndex+=(curLen);
			
			#add end
			cPIndex.upto(curProfile.length()-1){|i|
				fullProfile.push(curProfile[i].to_s());
			}
		else
			cPIndex.upto(curProfile.length()-1){|i|
				fullProfile.push(curProfile[i].to_s());
			}
		end
		
		
		#now add NaNs to end if appropriate
		fEnd.upto(idealEnd-1-strandOffsetSt){|i| #NaNs at end #updated 20130130
			fullProfile.push(NaN);
		}
		if strand=="-"
			fullProfile.reverse!();
		end
		return fullProfile;
	end

	private :doThisGFF;

	def profile(outFP)
		outFile = smartGZOut(outFP, "w");
		outFile.print(self.header());
		@allGFFs.each{|chr, st, en, strand, name|
			if @chromosomes.hasKey(chr) #added 20130323 to avoid printing missing gff entries for missing chrs
				begin
					fullProfile = doThisGFF(chr, st, en, strand);
				rescue Exception => e
					warn("Error raised for chr=#{chr} st=#{(st+1).to_s()}, en=#{(en+1).to_s}, strand=#{strand}"); 
					raise(e.message+"\n"+e.backtrace.join("\n")+"\n");
				end
				outFile.print(name+"\t"+fullProfile.join("\t")+"\n");
			end
		}
		outFile.close();
	end
	
	def revProfile(outFP)
		outFile = smartGZOut(outFP, "w");
		outFile.print(self.header());
		@allGFFs.each{|chr, st, en, strand, name|
			if @chromosomes.hasKey(chr) #added 20130323 to avoid printing missing gff entries for missing chrs
				fullProfile = doThisGFF(chr, st, en, strand).reverse();
				outFile.print(name+"\t"+fullProfile.join("\t")+"\n");
			end
		}
		outFile.close();
	end

	def makeProfile(chr, st, en, strand) #TODO: is this zero indexed?
		theProfile = [];
		if strand=="+"
			st.upto(en){|curStart|
				theProfile.push(@myScoringFtn.score(@chromosomes[chr], curStart));
			}
		elsif strand=="-"
			st.upto(en){|curStart|
				theProfile.push(@myScoringFtn.rcScore(@chromosomes[chr], curStart));
			}
		else
			raise("BAAAAAAAAAAAAAAD STRAND: "+[strand, chr, st, en].join("\t"));
		end
		return theProfile;	
	end
	def header()
		theHeader=[];
		0.upto(@flank-1){|i|
			theHeader.push("Up."+i.to_s());
		}
		if @scaleCentres
			0.upto(@avgFeatureLength-1){|i|
				theHeader.push("Mid."+i.to_s());
			}
		else
			0.upto(@maxFeatureLength-1){|i|
				theHeader.push("Mid."+i.to_s());
			}
		end
		0.upto(@flank-1){|i|
			theHeader.push("Down."+i.to_s());
		}
		return "\t"+theHeader.join("\t")+"\n";
	end

end

class GFFProfilerTrack < GFFProfiler
	def initialize(gffFP, genomeFP, flankingSeq, tracks)
		require "GENOMEDATA.rb";
		super(gffFP, genomeFP, "1", flankingSeq); #placeholder for scoring ftn to have length = 1
		#initialize(gffFP, genomeFP, myScoringFtn, flankingSeq)
		tracks = tracks.split(":");
		@trackP = loadMyFile(tracks[0]);
		@trackN=nil;
		if tracks.length()==2
			@trackN = loadMyFile(tracks[1]);
		end
	end

	def makeProfile(chr, st, en, strand)
		theProfile = [];
		if @trackP[chr].nil? || ((!@trackN.nil?) &&@trackN[chr].nil?)
			st.upto(en){|curStart|
				theProfile.push(NaN);
			}
			return theProfile;
		end
			
		begin
			if strand=="+" || @trackN.nil?
				st.upto(en){|curStart|
					theProfile.push(@trackP[chr][curStart]); #TODO: subtract one from coord?
				}
			elsif strand=="-"
				st.upto(en){|curStart|
					theProfile.push(@trackN[chr][curStart]);
				}
			else
				raise("BAAAAAAAAAAAAAAD STRAND: "+[strand, chr, st, en].join("\t"));
			end
		rescue NoMethodError =>e
			warn(e.message);
			warn(e.backtrace.inspect);
			raise ("Tried to access missing chr: #{chr}, #{st}, #{en}, #{strand}, [#{@trackP.keys().join(", ")}]");
		end
		theProfile.map!(){|e| 
			NaN if e.nil?; 
			e;
		}
		return theProfile;	
	end
end

class GFFProfilerWIG < GFFProfilerTrack
	def initialize(gffFP, genomeFP, flankingSeq, tracks)
		super(gffFP, genomeFP, flankingSeq, tracks);
	end
	def loadMyFile(file)
		#curData = makeNewAllData(@chromosomes, NaN);#commented out because then it takes forever and rusn out of memory.
		curData = readWigFile(curData,file);
		warn("Input all data");
		return curData;
	end
end

class GFFProfilerBedGraph < GFFProfilerTrack
	def initialize(gffFP, genomeFP, flankingSeq, tracks)
		super(gffFP, genomeFP, flankingSeq, tracks);
	end
	def loadMyFile(file)
		curData = makeNewAllData(@chromosomes, NaN);
		curData = readBedGraphFile(curData,file);
		return curData;
	end
end

class RangeScorerGFFProfiler <GFFProfiler
	def initialize(gffFP, genomeFP, myScoringFtn, flankingSeq, segLen)
		@myScorer = myScoringFtn;
		@segLen = segLen;
		super(gffFP, genomeFP, "1", flankingSeq); #placeholder for scoring ftn to have length = 1
	end
	
	def makeProfile(chr, st, en, strand)
		theProfile = [];
		sectionLen = (1.0+en-st);
		numSegs = (sectionLen/@segLen).round();
		1.upto(numSegs){|segNum|
			curStart =  ((sectionLen*(segNum-1))/numSegs).round()+st;
			curEnd =  (((sectionLen*segNum)/numSegs).round()+st-1);
			if curEnd>en
				raise("going past end: "+[chr,st,en,strand,curEnd,segNum,@segLen, numSegs, sectionLen].join(" ")+"\n")
			end
			curScore = @myScorer.scoreThisSeq(@chromosomes[chr], curStart, curEnd, strand);
			curScore = curScore[0].to_f();
			curStart.upto(curEnd){|i|
				theProfile.push(curScore);
			}
		}
		return theProfile;	
	end
end

class NoScaleGFFProfiler < GFFProfiler
	def initialize(gffFP, genomeFP, myScoringFtn, flankingSeq)
		super(gffFP, genomeFP, myScoringFtn, flankingSeq);
		maxLength=0;
		File.foreach(gffFP){|line|
			line.chomp!();
			if line==nil||line==""||line[0,1]=="#"||line[0,1]==">"||line[0,11]=="track name="
				next;
			end
			allData = line.split("\t");
			#chr start end strand name
			st = allData[3].to_i()-1;
			en = allData[4].to_i()-1;
			fLen=en-st+1;
			maxLength = [maxLength, fLen].max()
		}
		@avgFeatureLength = maxLength;
	end
	def doThisGFF(chr, st, en, strand)
		#print([chr, st, en, strand].join(" "));
		fStart = [0,st-@flank].max();
		fEnd = [@chromosomes[chr].length()-1, en+@flank].min();
		fEnd-= (@myScoringFtn.length()-1);
		#print(" "+[fStart, fEnd, @chromosomes[chr].length].join(" ")+"\n");
		curProfile = self.makeProfile(chr, fStart, fEnd, strand);
		fullProfile = [];
		#now add NaNs to front if appropriate
		strandOffset=0;
		if strand=="-"
			strandOffset = @myScoringFtn.length()-1;
		end
		(st-@flank).upto(fStart-1+strandOffset){|i| #NaNs at start
			fullProfile.push(NaN);
		}
		cPIndex=0;
		(fStart).upto(st-1-strandOffset){|i| #non-NaN start
			fullProfile.push(curProfile[cPIndex]);
			cPIndex+=1;
		}
		#middle
		curLen=[(en - st+1),(fEnd-st+1)].min();
		cPIndex.upto(cPIndex+curLen-1){|i|
			fullProfile.push(curProfile[i]);
		}
		cPIndex+=(curLen);
		
		#add end
		cPIndex.upto(curProfile.length()-1){|i|
			fullProfile.push(curProfile[i].to_s());
		}
		
		#now add NaNs to end if appropriate
		fEnd.upto(en+@flank-1-strandOffset){|i| #NaNs at end
			fullProfile.push(NaN);
		}
		if strand=="-"
			fullProfile.reverse!();
		end
		return fullProfile;
	end
end

class GFFProfilerEFBS < GFFProfiler
	def initialize(gffFP, genomeFP, myScoringFtn, flankingSeq, foldedFP)
		super(gffFP, genomeFP, myScoringFtn, flankingSeq);
		@acc =Hash.new();
		File.foreach(foldedFP){|line|
			chr, path = line.chomp.split("\t");
			@acc[chr+"+"] = [0.0]*@chromosomes[chr].length();
			i=0;
			File.foreach(path+".f.plout"){|sline|
				sline.chomp!();
				if sline==nil||sline==""||sline[0,1]=="#"
					next;
				end
				pos, probUnpaired = sline.split(" ");
				@acc[chr+"+"][i]=probUnpaired.to_f();
				i+=1;
			}
			@acc[chr+"-"] = [0.0]*@chromosomes[chr].length();
			i=0;
			File.foreach(path+".r.plout"){|sline|
				sline.chomp!();
				if sline==nil||sline==""||sline[0,1]=="#"
					next;
				end
				pos, probUnpaired = sline.split(" ");
				@acc[chr+"-"][i] = probUnpaired.to_f();
				i+=1;
			}
		}
		print ("Done initializing\n");
	end
	def makeProfile(chr, st, en, strand)
		curProfile = super(chr, st, en, strand);
		curProbUnpaired=1.0;
		index=0;
		myPArray = @acc[chr+strand];
		0.upto(curProfile.length()-1){|j|#st.upto(en-@myScoringFtn.length()+1){|j|
			curPU = 1.0;
			0.upto(@myScoringFtn.length()-1){|k|
				if st+j+k>myPArray.length-1
					raise("ERROR: "+[chr, st, en, strand, st, j, k, st+j+k, myPArray.length, curPU, curProfile.length()].join(" ")+"\n");
				end
				curPU=curPU*myPArray[st+j+k];
			}
			if @myScoringFtn.is_a?(PWMScoringFtn)
				curProfile[j] = 0.0+(curProfile[j]+(Math.log(curPU)/Math.log(2)));
			else
				curProfile[j] = 0.0+(curProfile[j]*curPU);
			end
		}
		
		return curProfile;
	end
end



def getProfiler(gffFP, details)
	stuff = details.split("_");
	genome=stuff[0];
	accessibility=stuff[1];
	scoreMethod = stuff[2];
	motif = stuff[3, stuff.length()-1-3].join("_");
	flankingSeq=stuff[stuff.length()-1].to_i();
	#make Scoring method
	useAvg=false;
	if scoreMethod[-2,2]=="-A"
		scoreMethod = scoreMethod[0,scoreMethod.length()-2];
		useAvg=true;
	end
	case scoreMethod
	when "C"
		myScoringFtn = ConsensusScoringFtn.new(motif)
	when "GSF"
		require "GOMERSCORER.rb";
		myScoringFtn = GomerScorerFromFile.new(motif, 0);
	when "GST"
		require "GOMERSCORER.rb";
		myScoringFtn = GomerScorerFromFile.new(motif, 1);
	when "GSR"
		require "GOMERSCORER.rb";
		myScoringFtn = GomerScorerFromFile.new(motif, -1);
	when "D"
		myScoringFtn = NonSSConsensusScoringFtn.new(motif)
	when "NSSPLPMF"
		myScoringFtn = NonSSPLPMFromFileScoringFtn.new(motif)
	when "PLPMF"
		myScoringFtn = PLPMFromFileScoringFtn.new(motif)
	when "P"
		myScoringFtn = PWMFromIUPACScoringFtn.new(motif)
	when "N"
		myScoringFtn = NeutralPWMFromFileScoringFtn.new(motif);
	when "F"
		myScoringFtn = PWMFromFileScoringFtn.new(motif);
	when "Q"
		myScoringFtn = NonSSPWMFromIUPACScoringFtn.new(motif)
	when "G"
		myScoringFtn = NonSSPWMFromFileScoringFtn.new(motif);
	when "T"
		myScoringFtn = TestScoringFtn.new();
	when "B"
		myScoringFtn = ONLYBESTPWMFromFileScoringFtn.new(motif);
	end
	#get genome
	case genome
	when "MMv37-chr1"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr1.map";
	when "MMv37-chr2"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr2.map";
	when "MMv37-chr3"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr3.map";
	when "MMv37-chr4"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr4.map";
	when "MMv37-chr5"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr5.map";
	when "MMv37-chr6"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr6.map";
	when "MMv37-chr7"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr7.map";
	when "MMv37-chr8"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr8.map";
	when "MMv37-chr9"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr9.map";
	when "MMv37-chr10"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr10.map";
	when "MMv37-chr11"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr11.map";
	when "MMv37-chr12"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr12.map";
	when "MMv37-chr13"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr13.map";
	when "MMv37-chr14"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr14.map";
	when "MMv37-chr15"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr15.map";
	when "MMv37-chr16"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr16.map";
	when "MMv37-chr17"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr17.map";
	when "MMv37-chr18"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr18.map";
	when "MMv37-chr19"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchr19.map";
	when "MMv37-chrX"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchrX.map";
	when "MMv37-chrY"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/justchrY.map";
	when "MMv37"
		genomeFP=HOME+"/genomes/mouse/NCBIv37/chr.map";
	when "Dmel5NH"
		genomeFP=HOME+"/genomes/drosophila/Dmel_Release5/chr_noHet.map";
	when "Dmel5"
		genomeFP=HOME+"/genomes/drosophila/Dmel_Release5/chr.map";
	when "hg18-chr1"
		genomeFP=HOME+"/genomes/human/hg18/justchr1.map";
	when "hg18-chr2"
		genomeFP=HOME+"/genomes/human/hg18/justchr2.map";
	when "hg18-chr3"
		genomeFP=HOME+"/genomes/human/hg18/justchr3.map";
	when "hg18-chr4"
		genomeFP=HOME+"/genomes/human/hg18/justchr4.map";
	when "hg18-chr5"
		genomeFP=HOME+"/genomes/human/hg18/justchr5.map";
	when "hg18-chr6"
		genomeFP=HOME+"/genomes/human/hg18/justchr6.map";
	when "hg18-chr7"
		genomeFP=HOME+"/genomes/human/hg18/justchr7.map";
	when "hg18-chr8"
		genomeFP=HOME+"/genomes/human/hg18/justchr8.map";
	when "hg18-chr9"
		genomeFP=HOME+"/genomes/human/hg18/justchr9.map";
	when "hg18-chr10"
		genomeFP=HOME+"/genomes/human/hg18/justchr10.map";
	when "hg18-chr11"
		genomeFP=HOME+"/genomes/human/hg18/justchr11.map";
	when "hg18-chr12"
		genomeFP=HOME+"/genomes/human/hg18/justchr12.map";
	when "hg18-chr13"
		genomeFP=HOME+"/genomes/human/hg18/justchr13.map";
	when "hg18-chr14"
		genomeFP=HOME+"/genomes/human/hg18/justchr14.map";
	when "hg18-chr15"
		genomeFP=HOME+"/genomes/human/hg18/justchr15.map";
	when "hg18-chr16"
		genomeFP=HOME+"/genomes/human/hg18/justchr16.map";
	when "hg18-chr17"
		genomeFP=HOME+"/genomes/human/hg18/justchr17.map";
	when "hg18-chr18"
		genomeFP=HOME+"/genomes/human/hg18/justchr18.map";
	when "hg18-chr19"
		genomeFP=HOME+"/genomes/human/hg18/justchr19.map";
	when "hg18-chr20"
		genomeFP=HOME+"/genomes/human/hg18/justchr20.map";
	when "hg18-chr21"
		genomeFP=HOME+"/genomes/human/hg18/justchr21.map";
	when "hg18-chr22"
		genomeFP=HOME+"/genomes/human/hg18/justchr22.map";
	when "hg18-chrX"
		genomeFP=HOME+"/genomes/human/hg18/justchrX.map";
	when "hg18-chrY"
		genomeFP=HOME+"/genomes/human/hg18/justchrY.map";
	when "hg18"
		genomeFP=HOME+"/genomes/human/hg18/chr.map";
	when "test"
		genomeFP=HOME+"/genomes/test/chrX.map";
	when "20110203"
		genomeFP=HOME+"/genomes/sc/20110203_R64/chrX.map";
	when "20110203F"
		genomeFP=HOME+"/genomes/sc/20110203_R64/chr1-8.map";
	when "20110203L"
		genomeFP=HOME+"/genomes/sc/20110203_R64/chr9-16.map";
	when "2008T"
		genomeFP=HOME+"/genomes/sc/May2008/chrX.mapT";
		foldedFP = HOME+"/genomes/sc/May2008/folded.mapT";
	when "2008F"
		genomeFP=HOME+"/genomes/sc/May2008/chr1-8.map";
	when "2008L"
		genomeFP=HOME+"/genomes/sc/May2008/chr9-16.map";
	when "2008"
		genomeFP=HOME+"/genomes/sc/May2008/chrX.map";
		foldedFP = HOME+"/genomes/sc/May2008/folded.map";
	when "2007"
		genomeFP=HOME+"/genomes/sc/20070717/chrX.map";
	when "CA-A21"
		genomeFP=HOME+"/genomes/calb/20100624_A21/chrX.map";
	when "SP"
		genomeFP=HOME+"/genomes/S_pombe/chr.map";
	when "Ashbya-gossypii" 
		genomeFP=HOME+"/genomes/Ashbya_gossypii/chr.map"
	when "Aspergillus-nidulans" 
		genomeFP=HOME+"/genomes/Aspergillus_nidulans/chr.map"
	when "Candida-glabrata" 
		genomeFP=HOME+"/genomes/Candida_glabrata/chr.map"
	when "Debaryomyces-hansenii" 
		genomeFP=HOME+"/genomes/Debaryomyces_hansenii/chr.map";
	when "Kluyveromyces-lactis" 
		genomeFP=HOME+"/genomes/Kluyveromyces_lactis/chr.map";
	when "Kluyveromyces-thermotolerans" 
		genomeFP=HOME+"/genomes/Kluyveromyces_thermotolerans/chr.map";
	when "Neurospora-crassa" 
		genomeFP=HOME+"/genomes/Neurospora_crassa/chr.map";
	when "S-bayanus" 
		genomeFP=HOME+"/genomes/S_bayanus/chrX.map";
	when "S-castellii" 
		genomeFP=HOME+"/genomes/S_castellii/chrX.map";
	when "S-cerevisiae" 
		genomeFP=HOME+"/genomes/sc/May2008/chrX.map";
	when "S-mikatae-WashU" 
		genomeFP=HOME+"/genomes/S_mikatae_WashU/chrX.map";
	when "S-paradoxus" 
		genomeFP=HOME+"/genomes/S_paradoxus/chrX.map";
	when "Yarrowia-lipolytica" 
		genomeFP=HOME+"/genomes/Yarrowia_lipolytica/chr.map";
	when "Zygosaccharomyces-rouxii" 
		genomeFP=HOME+"/genomes/Zygosaccharomyces_rouxii/chr.map"
	else
		raise("Unrecognized genome!");
	end
	#make profiler
	if scoreMethod=="WIG"
		profiler = GFFProfilerWIG.new(gffFP, genomeFP, flankingSeq, motif); 
	elsif scoreMethod=="BEDGR"
		profiler = GFFProfilerBedGraph.new(gffFP, genomeFP, flankingSeq, motif); 
	elsif scoreMethod=="GSF" || scoreMethod=="GST"
		profiler = RangeScorerGFFProfiler.new(gffFP, genomeFP, myScoringFtn, flankingSeq, 20); 
	elsif accessibility=="Y"
		profiler = GFFProfilerEFBS.new(gffFP, genomeFP, myScoringFtn, flankingSeq, foldedFP);
	elsif accessibility=="N"
		profiler = GFFProfiler.new(gffFP, genomeFP, myScoringFtn, flankingSeq);
	else
		raise("BAD PARAMETER");
	end
	profiler.useAvg() if useAvg;
	return profiler;
end

def getScoringFunction(scoringFtn, fileOrMotif)
	theScoringFtn = case scoringFtn
		when "LookupScoringFtn" then LookupScoringFtn.new(fileOrMotif);
		when "GCContentScoringFtn" then GCContentScoringFtn.new()
		when "GContentScoringFtn, CContentScoringFtn" then [GContentScoringFtn.new(), CContentScoringFtn.new()];
		when "AContentScoringFtn, TContentScoringFtn" then [AContentScoringFtn.new(), TContentScoringFtn.new()];
		when "AContentScoringFtn" then AContentScoringFtn.new();
		when "TContentScoringFtn" then TContentScoringFtn.new();
		when "GContentScoringFtn" then GContentScoringFtn.new();
		when "CContentScoringFtn" then CContentScoringFtn.new();
		when "PWMFromFileScoringFtn" then PWMFromFileScoringFtn.new(fileOrMotif);
		when "NonSSPWMFromFileScoringFtn" then NonSSPWMFromFileScoringFtn.new(fileOrMotif);
		when "NonSSConsensusScoringFtn" then NonSSConsensusScoringFtn.new(fileOrMotif);
		when "ConsensusScoringFtn" then ConsensusScoringFtn.new(fileOrMotif);
		when "" then fileOrMotif;
		else raise("Unknown scoring function: "+scoringFtn+"\n");
	end
	return theScoringFtn;
end

def getScorer_noFold(scorer, theScoringFtn)
	theScorer = case scorer
		when "BendScorer" then BendScorer.new(theScoringFtn);
		when "SimpleAverage" then SimpleAverage.new(theScoringFtn);
		when "MeltingTemp" then MeltingTemp.new();
		when "RatioScorer(SimpleSum)" then RatioScorer.new(SimpleSum.new(theScoringFtn[0]), SimpleSum.new(theScoringFtn[1]));
		when "SimpleSum" then SimpleSum.new(theScoringFtn);
		when "GomerScorerFromFile(<file>,true)" then GomerScorerFromFile.new(theScoringFtn,1);
		when "GomerScorerFromFile(<file>,false)" then GomerScorerFromFile.new(theScoringFtn,0);
		when "GomerScorerFromFile(<file>,1)" then GomerScorerFromFile.new(theScoringFtn,1);
		when "GomerScorerFromFile(<file>,-1)" then GomerScorerFromFile.new(theScoringFtn,-1);
		when "GomerScorerFromFile(<file>,0)" then GomerScorerFromFile.new(theScoringFtn,0);
		#when "AdditivePUnpairedScorer" then AdditivePUnpairedScorer.new(theScoringFtn, pUnPMap);
		when "HairpinLoopScorer" then HairpinLoopScorer.new(theScoringFtn);
		else raise("Unknown scorer: #{scorer} #{scoringFtn}\n");
	end
	return theScorer;
end

def getScorer(name, fileOrMotif, scorer, scoringFtn, chromosomes, pUnPMap)
	theScoringFtn = case scoringFtn
		when "LookupScoringFtn" then LookupScoringFtn.new(fileOrMotif);
		when "ReverseLookupScoringFtn" then ReverseLookupScoringFtn.new(fileOrMotif);
		when "GCContentScoringFtn" then GCContentScoringFtn.new()
		when "GContentScoringFtn, CContentScoringFtn" then [GContentScoringFtn.new(), CContentScoringFtn.new()];
		when "AContentScoringFtn, TContentScoringFtn" then [AContentScoringFtn.new(), TContentScoringFtn.new()];
		when "AContentScoringFtn" then AContentScoringFtn.new();
		when "TContentScoringFtn" then TContentScoringFtn.new();
		when "GContentScoringFtn" then GContentScoringFtn.new();
		when "CContentScoringFtn" then CContentScoringFtn.new();
		when "PWMFromFileScoringFtn" then PWMFromFileScoringFtn.new(fileOrMotif);
		when "NonSSConsensusScoringFtn" then NonSSConsensusScoringFtn.new(fileOrMotif);
		when "ConsensusScoringFtn" then ConsensusScoringFtn.new(fileOrMotif);
		when "" then nil;
		else raise("Unknown scoring function: #{name} #{fileOrMotif} #{scorer} #{scoringFtn}\n");
	end
	theScorer = case scorer
		when "BendScorer" then BendScorer.new(theScoringFtn);
		when "SimpleAverage" then SimpleAverage.new(theScoringFtn);
		when "SimpleAverageSym" then SimpleAverageSym.new(theScoringFtn);
		when "MeltingTemp" then MeltingTemp.new();
		when "PolyNScorer" then PolyNScorer.new(fileOrMotif);
		when "SSPolyNScorer" then SSPolyNScorer.new(fileOrMotif);
		when "GQuad" then GQuad.new();
		when "RatioScorer(SimpleSum)" then RatioScorer.new(SimpleSum.new(theScoringFtn[0]), SimpleSum.new(theScoringFtn[1]));
		when "SimpleSum" then SimpleSum.new(theScoringFtn);
		when "SimpleSumSym" then SimpleSumSym.new(theScoringFtn);
		when "GomerScorerFromFile(<file>,true)" then GomerScorerFromFile.new(fileOrMotif,1);
		when "GomerScorerFromFile(<file>,false)" then GomerScorerFromFile.new(fileOrMotif,0);
		when "GomerScorerFromFile(<file>,1)" then GomerScorerFromFile.new(fileOrMotif,1);
		when "GomerScorerFromFile(<file>,-1)" then GomerScorerFromFile.new(fileOrMotif,-1);
		when "GomerScorerFromFile(<file>,0)" then GomerScorerFromFile.new(fileOrMotif,0);
		when "BestPUnpairedScorer" then BestPUnpairedScorer.new(theScoringFtn, pUnPMap);
		when "AdditivePUnpairedScorer" then AdditivePUnpairedScorer.new(theScoringFtn, pUnPMap);
		when "HairpinLoopScorer" then HairpinLoopScorer.new(theScoringFtn);
		when "AverageWigScorer" then AverageWigScorer.new(fileOrMotif);
		when "SumWigScorer" then SumWigScorer.new(fileOrMotif);
		else raise("Unknown scorer: #{name} #{fileOrMotif} #{scorer} #{scoringFtn}\n");
	end
	return theScorer;
end

def getScorer_noPunPByFtn(name, fileOrMotif, scorer, scoringFtn)
	theScoringFtn = case scoringFtn
		when "LookupScoringFtn" then LookupScoringFtn.new(fileOrMotif);
		when "GCContentScoringFtn" then GCContentScoringFtn.new()
		when "GContentScoringFtn, CContentScoringFtn" then [GContentScoringFtn.new(), CContentScoringFtn.new()];
		when "AContentScoringFtn, TContentScoringFtn" then [AContentScoringFtn.new(), TContentScoringFtn.new()];
		when "AContentScoringFtn" then AContentScoringFtn.new();
		when "TContentScoringFtn" then TContentScoringFtn.new();
		when "GContentScoringFtn" then GContentScoringFtn.new();
		when "CContentScoringFtn" then CContentScoringFtn.new();
		when "PWMFromFileScoringFtn" then PWMFromFileScoringFtn.new(fileOrMotif);
		when "NonSSConsensusScoringFtn" then NonSSConsensusScoringFtn.new(fileOrMotif);
		when "ConsensusScoringFtn" then ConsensusScoringFtn.new(fileOrMotif);
		when "" then nil;
		else raise("Unknown scoring function: #{name} #{fileOrMotif} #{scorer} #{scoringFtn}\n");
	end
	theScorer = case scorer
		when "BendScorer" then BendScorer.new(theScoringFtn);
		when "SimpleAverage" then SimpleAverage.new(theScoringFtn);
		when "MeltingTemp" then MeltingTemp.new();
		when "RatioScorer(SimpleSum)" then RatioScorer.new(SimpleSum.new(theScoringFtn[0]), SimpleSum.new(theScoringFtn[1]));
		when "SimpleSum" then SimpleSum.new(theScoringFtn);
		when "GomerScorerFromFile(<file>,true)" then GomerScorerFromFile.new(fileOrMotif,1);
		when "GomerScorerFromFile(<file>,false)" then GomerScorerFromFile.new(fileOrMotif,0);
		when "GomerScorerFromFile(<file>,1)" then GomerScorerFromFile.new(fileOrMotif,1);
		when "GomerScorerFromFile(<file>,-1)" then GomerScorerFromFile.new(fileOrMotif,-1);
		when "GomerScorerFromFile(<file>,0)" then GomerScorerFromFile.new(fileOrMotif,0);
		when "HairpinLoopScorer" then HairpinLoopScorer.new(theScoringFtn);
		when "PolyNScorer" then PolyNScorer.new(fileOrMotif);
		else raise("Unknown scorer: #{name} #{fileOrMotif} #{scorer} #{scoringFtn}\n");
	end
	return theScorer;
end
