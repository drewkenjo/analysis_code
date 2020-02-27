package run

class Run {
	def run_number = 3000 // default run number
    def polarity = true //

    // characterize run numbers 
    def inbending_low_lumi=[4887, 4888, 4900, 4893, 4895, 5418, 5419]
    def outbending_low_lumi=[5443, 5444, 5453]

    def polarity_by_run = {run->
    	def polarity
    	if (inbending_low_lumi.contains(run)) polarity = true
    	if (outbending_low_lumi.contains(run)) polarity = false
    	return polarity
    }
}