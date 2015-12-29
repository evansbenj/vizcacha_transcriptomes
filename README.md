# vizcacha_transcriptomes

This is a repository for the vizcacha rat transcriptome project.  There are three goals of the first stage of this project:

* assemble and compare Tympanoctomys and Octotomys transcriptomes from liver, heart, muscle, gonad, lung, and kidney.  We have one biological replicate for each species so this will be done by summing the difference between orthologs for each tissue and then seeing which one is most diverged.  We expect gonad to be most diverged because we have one male and one female.  And then we expect kidney to be most diverged because Tympa is arid adapted.

* BLAST each transcriptome assembly with a query of mammalian repetitive elements.  This will test the hypothesis that Tympa has a large genome due to the expansion of repetitive elements

* Combine the data with the transcriptome of tuco tucos to identify orthologs and paralogs.  There are several ways to do this; one might be to pool the tuco tuco and tympa transcriptomes and then save the top hits (at least top two) with octomys as a query. The expectation if Tympa is polyploid is that the top two would be Tympa sequences.  If Tympa is diploid, we would expect almost all of the top two hits to include one Tympa and one tuco tuco.
