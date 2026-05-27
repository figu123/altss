import socketserver

from RangeHTTPServer import RangeRequestHandler


class CORSAndRangeHandler(RangeRequestHandler):
    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Range, Content-Type")
        self.send_header(
            "Access-Control-Expose-Headers", "Content-Range, Accept-Ranges"
        )
        super().end_headers()

    def do_OPTIONS(self):
        self.send_response(200)
        self.end_headers()


if __name__ == "__main__":
    porta = 8080
    print(f"Servidor IGV rodando na porta {porta}!")
    httpd = socketserver.TCPServer(("", porta), CORSAndRangeHandler)
    httpd.serve_forever()
