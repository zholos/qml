# Remove example program.

/START OF CONMAX SUBPROGRAMS/,$!d


# No need to pass dimensions since we pass a struct pointer in PTTBL.

s/PTTBL(IPTB,INDM)/PTTBL(*)/g
s/IPTB,//g;s/INDM,//g


# Allow setting additional parameters when calling these functions directly.

/^..... *SUBROUTINE MULLER(/{
    h;s/(.*/(LIMMUL,NSRCH,/p
    g;s/.*(/     */
}
/^..... *CALL MULLER(/{
    h;s/(.*/(5,NSRCH,/p
    g;s/.*(/     */
}
/^..... *LIMMUL=5 *$/d

/^..... *SUBROUTINE SEARSL(/{
    h;s/(.*/(INITLM,NADD,LIMS1,/p
    g;s/.*(/     */
}
/^..... *CALL SEARSL(/{
    h;s/(.*/(6,4,LIMS1,/p
    g;s/.*(/     */
}
/^..... *INITLM=6 *$/d
/^..... *NADD=4 *$/d
