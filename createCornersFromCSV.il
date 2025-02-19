; Get session and database handles
sess = axlGetWindowSession()
sdb = axlGetMainSetupDB(sess)
N=12
; Main function: Process CSV file containing multiple binary codes
procedure(createCornersFromCSV(csvFileName)
    let((inPort line cornerCount)
        ; Try to open the CSV file
        when((inPort = infile(csvFileName)) == nil
            error("Failed to open file %s" csvFileName)
            return(nil)
        )
        
        cornerCount = 0  ; Track number of corners created
        
        ; Read CSV file line by line
        while(gets(line inPort)
            ; Clean whitespace characters from the line
            line = stripString(line)
            
            ; Check if it's a valid N-bit binary code
            when(and(line != "" isValidBinary(line))
                ; Create new corner using binary code as name
                corn = axlPutCorner(sdb line)
                
                ; Set variables for this corner
                setBinaryToVars(line corn)
                
                cornerCount = cornerCount + 1
                printf("Successfully created corner: %s\n" line)
            )
        )
        
        ; Close file
        close(inPort)
        
        ; Return statistics
        printf("Total corners created: %d\n" cornerCount)
        cornerCount
    )
)

; Helper function: Validate N-bit binary code
procedure(isValidBinary(str)
    let((valid)
        valid = (strlen(str) == N)  ; Check length is N
        
        ; Verify only contains 0 and 1
        when(valid
            for(i 1 N
                unless(member(substring(str i 1) list("0" "1"))
                    valid = nil
                    break()
                )
            )
        )
        
        unless(valid
            printf("Warning: Invalid binary code '%s' - skipped\n" str)
        )
        valid
    )
)

; Set variables for a single corner
procedure(setBinaryToVars(binString corn)
    let((varNames)
        ; Define array of variable names
        varNames = list("VC1I0" "VC1I1" "VC1I2" "VC1I3" "VC1I4" "VC1I5" 
                       "VC1Q0" "VC1Q1" "VC1Q2" "VC1Q3" "VC1Q4" "VC1Q5" )
        
        ; Process each bit in binary string
        for(i 1 N
            ; Get current bit value ("0" or "1")
            binValue = substring(binString i 1)
            if(binValue == "1" then
                binValue = "1.2"
            else
                binValue = "0")
            ; Get corresponding variable name (nth is 0-based)
            varName = nth(i-1 varNames)
            
            ; Set the binary value to corresponding variable
            axlPutVar(corn varName binValue)
        )
    )
)

; Helper function: Clean whitespace from string
procedure(stripString(str)
    rexCompile("^[ \t\n\r]*")
    str = rexReplace(str "" 0)
    rexCompile("[ \t\n\r]*$")
    rexReplace(str "" 0)
)

; Usage example:
; createCornersFromCSV("path/to/your/corners.csv")