# State    

```{mermaid}

stateDiagram-v2
    
    1 : FILE
    2 : Just DATA
    3 : Just SAVE
    4 : LOOP
    5 : NAME
    6 : JUST SAVE LOOP
    7 : SAVE NAME
    8 : LOOP NAME
    9 : DATA
    10 : SAVE LOOP NAME
    11 : SAVE
    12 : LOOP VALUE
    13 : SAVE LOOP VALUE

    [*] --> 1
    
    1 --> 2 : DATA
    2 --> 3 : SAVE
    2 --> 4 : LOOP
    2 --> 5 : NAME
    3 --> 6 : LOOP
    3 --> 7 : NAME
    4 --> 8 : NAME
    5 --> 9 : VALUE
    6 --> 10 : NAME
    7 --> 11 : VALUE
    8 --> 8 : NAME
    8 --> 12 : VALUE
    9 --> 2 : DATA
    9 --> 3 : SAVE
    9 --> 4 : LOOP
    9 --> 5 : NAME
    10 --> 10 : NAME
    10 --> 13 : VALUE
    11 --> 9 : SAVE_END
    11 --> 6 : LOOP
    11 --> 7 : NAME
    12 --> 2 : DATA
    12 --> 3 : SAVE
    12 --> 4 : LOOP
    12 --> 5 : NAME
    12 --> 12 : VALUE
    13 --> 9 : SAVE_END
    13 --> 6 : LOOP
    13 --> 7 : NAME
    13 --> 13 : VALUE

    9 --> [*] : EOF
    11 --> [*] : EOF
    12 --> [*] : EOF
    13 --> [*] : EOF

```
State diagram of an mmCIF file parser. 