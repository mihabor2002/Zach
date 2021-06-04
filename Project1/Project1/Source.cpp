#include "Client.h"

int main(int argc, char** argv)
{
	Client c;

	const char* serverName;
	serverName = (argc < 2) ? "127.0.0.1" : argv[1];

	int err;
	int sock;
	struct sockaddr_in server_addr;
	//struct hostent *hostinfo;

	// èíèöèàëèçàöèÿ windows sockets
	WSADATA wsaData;
	if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) {
		printf("WSAStartup failed\n");
		return -1;
	}

	// Ïîëó÷àåì èíôîðìàöèþ î ñåðâåðå ïî åãî DNS èìåíè
	// èëè òî÷å÷íîé íîòàöèè IP àäðåñà.
	//hostinfo = gethostbyname(SERVER_NAME);
	//if (hostinfo == NULL) {
	// fprintf(stderr, "Unknown host %s.\n", SERVER_NAME);
	// exit(EXIT_FAILURE);
	//}
	// ìîæíî áûëî áû èñïîëüçîâàòü GetAddrInfo()

	// Çàïîëíÿåì àäðåñíóþ ñòðóêòóðó äëÿ ïîñëåäóþùåãî
	// èñïîëüçîâàíèÿ ïðè óñòàíîâëåíèè ñîåäèíåíèÿ
	server_addr.sin_family = PF_INET;
	server_addr.sin_port = htons(SERVER_PORT);
	//server_addr.sin_addr = *(struct in_addr*) hostinfo->h_addr;
	unsigned int iaddr;
	inet_pton(AF_INET, serverName, &iaddr);
	server_addr.sin_addr.s_addr = iaddr;

	// Ñîçäàåì TCP ñîêåò.
	sock = socket(PF_INET, SOCK_STREAM, 0);
	if (sock < 0) {
		perror("Client: socket was not created");
		exit(EXIT_FAILURE);
	}

	// Óñòàíàâëèâàåì ñîåäèíåíèå ñ ñåðâåðîì
	err = connect(sock, (struct sockaddr*)&server_addr, sizeof(server_addr));
	if (err < 0) {
		perror("Client: connect failure");
		exit(EXIT_FAILURE);
	}
	fprintf(stdout, "Connection is ready\n");

	// Îáìåíèâàåìñÿ äàííûìè
	while (1)
	{
		if (c.write_to_server(sock) < 0) break;
		if (c.read_from_server(sock) < 0) break;
	}
	fprintf(stdout, "The end\n");

	// Çàêðûâàåì socket
	closesocket(sock);
	WSACleanup();
	exit(EXIT_SUCCESS);
}