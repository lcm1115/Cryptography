def exp(a, b, data):
    p = data[0]
    mod = 1
    mul = a

    while b > 0:
        if b & 1:
            mod *= mul
            mod %= p
        mul *= mul
        mul %= p
        b >>= 1

    return mod

def ver(value, witness, element, data):
    return value == exp(witness, element, data)

def update(witness, element, data):
    return exp(witness, element, data)

p = 12836264062479034298175995644061129422150071440384970707002735673228505544198116212051300182327092825187408783724842480146482000632590656200300492111164887
q = 11681613303664159082395334899365418983487763493752910271941113595633508602924057867263172050641028359785090889903808808552693945098092117676177016998552123
n = 149948273041601231577286268576025719223767157169380463874000007363838290086835743831836498738961107020142153660800208582306590739133763852581936960971477296449288882798156128802789076738213965495324063300228274794684068279367772855739021658250697773479940019925577318702662220797005899727962211185641916905101
e = 65537
d = 83715418196120473353698327371228543275987862027106242163751564298571481697166990113417850258019865797335570140455297493294106360592722974840479293604916030680354852722712849154536159350298780377888364540963504705528090634016689991385136265262527255798579251877052506451160501435722882115360727201975544840397
g = 3 
