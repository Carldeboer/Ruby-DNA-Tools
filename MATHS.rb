
def smooth(values, smoothBy)
	window = [];
	windex = 0;
	wcount=0;
	total =0.0;
	newValues = [];
	before = ((-1.0+smoothBy)/2.0).ceil();
	after = ((-1.0+smoothBy)/2.0).floor();
	#p([before,after,smoothBy]);
	#preload
	0.upto(after-1){|i|
		wcount+=1;
		total = total+values[i];
	}
	#start
	0.upto(before){|i|
		total = total+values[i+after];
		wcount+=1;
		newValues[i] = total/wcount;
	}
	#p(wcount);
	#mid
	(before+1).upto(values.length()-after-1){|i|
		total+=values[i+after]-values[i-before-1];
		newValues[i]=total/smoothBy;
	}
	#end
	wcount = smoothBy;
	(values.length()-after).upto(values.length()-1){|i|
		total=total-values[i-before-1];
		wcount-=1;
		newValues[i]=total/wcount;
	}
	if newValues.length()!= values.length()
		raise("smoothing resulted in different sized array: "+ newValues.length().to_s()+" "+values.length().to_s()+"\n");
	end
	return newValues
end

VERBOSE=false;
def findPeaksDPWithMax(data,direction,penaltyBP) #returns an array of peaks where each peak =[start,score,end,max,maxpos];
	peaks = findPeaksDP(data, direction, penaltyBP);
	0.upto(peaks.length()-1){|i|
		theMax=-99999;
		theMaxPosS=nil;#start
		theMaxPosL=nil;#last
		peaks[i][0].upto(peaks[i][2]){|j|
			if (data[j]>theMax)
				theMax=data[j];
				theMaxPosS=j;
				theMaxPosL=j;
			elsif data[j]==theMax
				theMaxPosL=j;
			end
		}
		peaks[i].push(theMax);
		peaks[i].push((theMaxPosS+theMaxPosL)/2);
	}
	return peaks;
end

def findPeaksDP(data, direction, penaltyBP) #bppenalty should be positive
	score=[0.0]*data.length();
	origin=[0]*data.length();
	if direction>0
		i = 1;
		score[0]=data[0]-penaltyBP;
		origin[0]=0
	else
		i=data.length()-2;
		score[data.length()-1]=data[data.length()-1]-penaltyBP;
		origin[data.length()-1]=data.length()-1;
	end
	#p(score);
	#p(origin);
	#raise("BLAH");
	numBases=1;
	while (i>=0 && i<=data.length()-1)
		curScore=data[i]-penaltyBP;
		cumScore=score[i-direction]+curScore;
		if curScore>cumScore
			score[i]=curScore;
			origin[i]=i;
		else
			score[i]=cumScore;
			origin[i]=origin[i-direction];
		end
		#now find i for which score is max.
		i+=direction;
		numBases+=1;
	end
	if direction<0 #go in reverse order to find peaks
		i = 0;
	else
		i=data.length()-1;
	end
	allPeaks=[];
	while (i>=0 && i<=data.length()-1)
		if (score[i]>0)
			st= [origin[i],i].min();
			en= [origin[i],i].max();
			theMax=i;
			st.upto(en){|j|
				if score[j]>score[theMax];
					theMax=j;
				end
			}
			allPeaks.push([[origin[i],theMax].min(),score[theMax],[origin[i],theMax].max()]);
			i=origin[i];
		end
		i-=direction; #reversed from before
	end
	return allPeaks;
end

def findPeakDP(data, direction, penaltyBP,penaltyDist,degRate) #penalties should be positive (unless you want distant things, then penaltyDist can be negative)
	score=[0.0]*data.length();
	origin=[0]*data.length();
	if direction>0
		i = 1;
		score[0]=data[0]-penaltyBP;
		origin[0]=0
	else
		i=data.length()-2;
		score[data.length()-1]=data[data.length()-1]-penaltyBP;
		origin[data.length()-1]=data.length()-1;
	end
	#p(score);
	#p(origin);
	#raise("BLAH");
	maxI=0;
	maxScore=-999999.0;
	curDP = penaltyDist;
	numBases=1;
	while (i>=0 && i<=data.length()-1)
		curScore=data[i]-penaltyBP-(numBases*penaltyDist)**degRate;
		cumScore=score[i-direction]+curScore;
		if curScore>cumScore
			score[i]=curScore;
			origin[i]=i;
		else
			score[i]=cumScore;
			origin[i]=origin[i-direction];
		end
		#now find i for which score is max.
		if score[i]>maxScore
			maxScore=score[i];
			maxI=i;
		end
		i+=direction;
		numBases+=1;
	end
	return [[maxI,origin[maxI]].min(), maxScore,[maxI,origin[maxI]].max()];
end

def findPeakBoundaries(data, direction, theMin,theMax)
	found=false;
	
	if direction>0
		curPos = 0;
	else
		curPos=data.length()-1;
	end
	foundTop=false;
	while !foundTop && curPos>=0 && curPos<=data.length()-1
		if data[curPos]>theMax
			foundTop=true;
		end
		curPos+=direction;
	end
	if !foundTop
		print("Found no top\n") if VERBOSE;
		return [nil,nil,nil];
	end
	curPos=curPos-direction;
	theTop=curPos;
	print("Found top at #{theTop}\n") if VERBOSE;
	#find bottoms
	foundBottom=false;
	while !foundBottom && curPos>=0
		#if data[curPos].nil?
		#	p([curPos,theTop]);
		#end
		if data[curPos]>=theMin
			peakStart=curPos;
		else
			foundBottom=true;
		end
		curPos=curPos-1;
	end
	print("Found bottom at #{peakStart}\n") if VERBOSE;
	foundBottom=false;
	curPos=theTop;
	while !foundBottom && curPos<=data.length()-1
		if data[curPos]>=theMin
			peakEnd=curPos;
		else
			foundBottom=true;
		end
		curPos=curPos+1;
	end
	print("Found bottom at #{peakEnd}\n") if VERBOSE;
	return [peakStart,theTop,peakEnd];
end


def findPeak(data, start, dirs)
	found=false;
	
	peak=start; #rear (front)  |*    ||
	last = data[peak];
	di=0;
	while !found #look for peak
		if di>=dirs.length()
			break;
		elsif peak+dirs[di]>=data.length || peak+dirs[di]<0;
			di+=1;
		else
			this = data[peak+dirs[di]];
			if this>=last
				peak=peak+dirs[di];
				last=this;
			else
				di=di+1;
			end
		end
	end
	return peak;
end

def findTraugh(data, start, dirs)
	found=false;
	traugh=start; #rear (front)  |*    ||
	last = data[traugh];
	di=0;
	while !found #look for traugh
		if di>=dirs.length()
			break;
		elsif traugh+dirs[di]>=data.length || traugh+dirs[di]<0;
			di+=1;
		else
			this = data[traugh+dirs[di]];
			if this<=last
				traugh=traugh+dirs[di];
				last=this;
			else
				di=di+1;
			end
		end
	end
	return traugh;
end
