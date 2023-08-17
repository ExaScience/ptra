FROM golang:1.17 AS build-stage

WORKDIR /app

COPY go.mod go.sum ./
RUN go mod download

COPY . ./


RUN CGO_ENABLED=0 GOOS=linux go build -o /ptra

# Run the tests in the container
FROM build-stage AS run-test-stage
RUN go test -v ./...

# Deploy the application binary into a lean image
FROM gcr.io/distroless/base-debian11 AS build-release-stage

WORKDIR /

COPY --from=build-stage /ptra /ptra

USER nonroot:nonroot

ENTRYPOINT ["/ptra"]
