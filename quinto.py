"https://www.youtube.com/watch?v=FuEG31oLGBk"

def Inercia(b,d):
    Ixx=(pow(b,3)*d)/12;
    Iyy=(pow(d,3)*b)/12;
    Izz=Ixx+Iyy;
    print("Área Momento de inercia Ixx:", Ixx)
    print("Área Momento de inercia Iyy:", Iyy)
    print("Área Momento de inercia Izz:", Izz)
Area_momento_de_Inercia=Inercia(10,2)

