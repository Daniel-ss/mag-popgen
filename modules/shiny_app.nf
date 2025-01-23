process SHINY_APP {
    publishDir "${params.outdir}/shiny_app", mode: 'copy'

    afterScript 'echo "Workflow completed. Shiny app is now running."'

    input:
    path pogenom_results
    path metadata

    output:
    path "app/*"

    script:
    """
    # Create app directory
    mkdir -p app

    cp $metadata metadata.csv
    cp ${workflow.projectDir}/bin/shiny_app.R app/app.R
    
    # Copy POGENOM result files from all MAGs
    for dir in $pogenom_results
    do
        cp \$dir/*.fst.txt app/
        cp \$dir/*.intradiv.txt app/
    done

    # Create a startup script
    echo '#!/usr/bin/env Rscript
    library(shiny)
    runApp(appDir=".", port=3838, host="0.0.0.0")' > app/start.R

    #Rscript -e "shiny::runApp('app/app.R', port = 8787, host = '0.0.0.0')"

    # Run the Shiny app in the background
    #Rscript -e "shiny::runApp('shiny_app.R', port = 3838, host = '0.0.0.0')" &> shiny_app.log
    """
}
