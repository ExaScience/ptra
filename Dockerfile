FROM golang:1.17 AS build-stage

WORKDIR /app

COPY go.mod go.sum ./
RUN go mod download

COPY . ./
COPY .docker/entrypoint.sh start.sh

RUN CGO_ENABLED=0 GOOS=linux go build -o /ptra

# Run the tests in the container
FROM build-stage AS run-test-stage
RUN go test -v ./...

# Deploy the application binary into a lean image
FROM debian AS build-release-stage

# Install mcl package for clustering
# ref. https://debian.pkgs.org/11/debian-main-amd64/mcl_14-137+ds-9+b1_amd64.deb.html
RUN apt-get update && apt-get install -y mcl

WORKDIR /
RUN mkdir -p /input /output

COPY .docker/entrypoint.sh .
COPY --from=build-stage /ptra /ptra
RUN chmod +x entrypoint.sh

CMD ["./entrypoint.sh"]
