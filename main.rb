"""
Energia cinética inicial ~1MeV. Iões lançados forward.
Depois de cada colisão há isotropia (temos de considerar conservação de momento linear e de energia – recuo do átomo com que colide, que inicialmente está em repouso)
Densidade do meio – N átomos /cm3 (number density à pressão atmosférica) PTN
A simulação para quando os iões têm energia final de ~1keV. Regista-se a posição final do ião.
Iões de Ar+ em atmosfera de Ar.
Usar como secção eficaz o modelo da esfera rígida (para determinar o caminho livre médio do ião)

Upgrades ao problema (pontos a adicionar resolvida a primeira aproximação):
1)	Os átomos do meio têm distribuição de velocidade normal com média em 3/2kT.
2)	Há um campo eléctrico uniforme aplicado.
3)	Electrões em vez de iões a serem lançados e uma secção eficaz mais realista que terá de ser tratada pelo métodos numéricos aprendidos.
"""

EcI = 1E9

Ecf = 1E3