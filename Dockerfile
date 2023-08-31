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
FROM alpine AS build-release-stage

WORKDIR /
RUN mkdir -p /input /output

COPY .docker/entrypoint.sh .
COPY --from=build-stage /ptra /ptra
RUN chmod +x entrypoint.sh

CMD ["./entrypoint.sh"]
